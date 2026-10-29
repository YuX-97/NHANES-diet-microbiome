# ============================================================
# dbRDA of dietary intake and the oral-rinse microbiome
# Overall population and oral conditions (Oral_CDC = 1, 2, 3, 4)
# Both sequential and marginal permutation tests are reported.
# Aitchison distance = Euclidean distance on CLR-transformed taxa.
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(openxlsx)
})

# 1. Input, output, and variables ------------------------------

input_file <- "dbRDA_analysis_data.tsv"
output_file <- "dbRDA_overall_oral_conditions_dual_tests.xlsx"

id_var <- "SEQN"
group_var <- "Oral_CDC"
energy_var <- "Energy_intake_kcal"

# Supply variable NAMES from the analysis dataset; no column positions are used.
# These two vectors are placeholders in this shared analysis script.
diet_vars <- character(0)      # Dietary intake variables (excluding energy)
clr_taxa_vars <- character(0)  # CLR-transformed genus variables

n_perm <- 1001
random_seed <- 3000
group_codes <- c("Overall", "1", "2", "3", "4")
group_labels <- c(
  Overall = "Overall",
  `1` = "Healthy",
  `2` = "Caries",
  `3` = "Periodontitis",
  `4` = "Coexistence"
)

# 2. Import data and apply energy adjustment ------------------

dat <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

stopifnot(length(diet_vars) > 0, length(clr_taxa_vars) > 0)
required_vars <- c(id_var, group_var, energy_var, diet_vars, clr_taxa_vars)
stopifnot(all(required_vars %in% names(dat)))

diet_raw <- as.data.frame(lapply(
  dat[, diet_vars, drop = FALSE],
  function(x) suppressWarnings(as.numeric(x))
))
energy <- suppressWarnings(as.numeric(dat[[energy_var]]))

# Restrict the overall sample to the four prespecified oral conditions.
# The analysis also requires complete dietary data and positive energy intake.
keep <- as.character(dat[[group_var]]) %in% c("1", "2", "3", "4") &
  complete.cases(diet_raw) & is.finite(energy) & energy > 0
df <- dat[keep, , drop = FALSE]
diet_raw <- diet_raw[keep, , drop = FALSE]
energy <- energy[keep]
stopifnot(nrow(df) > 0)

# Express food intake per 1,000 kcal, then z-standardize globally.
diet_ea <- as.data.frame(sweep(
  as.matrix(diet_raw), 1, energy / 1000, FUN = "/"
))
names(diet_ea) <- paste0("DIET_", seq_along(diet_vars))

# Drop globally invariant dietary variables before standardization.
global_sd <- vapply(diet_ea, sd, numeric(1), na.rm = TRUE)
global_valid <- is.finite(global_sd) & global_sd > 0
diet_ea <- diet_ea[, global_valid, drop = FALSE]
stopifnot(ncol(diet_ea) > 0)

diet_z <- as.data.frame(scale(diet_ea, center = TRUE, scale = TRUE))
diet_map <- data.frame(
  Internal = paste0("DIET_", seq_along(diet_vars)),
  Dietary_Variable = diet_vars,
  stringsAsFactors = FALSE
)
diet_map <- diet_map[diet_map$Internal %in% names(diet_z), , drop = FALSE]

# The upstream CLR preprocessing is assumed complete.
clr_all <- as.matrix(df[, clr_taxa_vars, drop = FALSE])
storage.mode(clr_all) <- "numeric"
stopifnot(!anyNA(clr_all), all(is.finite(clr_all)))
rownames(clr_all) <- as.character(df[[id_var]])

# 3. Helpers for model fitting and result extraction -----------

drop_linear_dependencies <- function(x) {
  x <- as.data.frame(x)
  if (ncol(x) <= 1) return(x)
  qr_fit <- qr(as.matrix(x))
  valid_columns <- sort(qr_fit$pivot[seq_len(qr_fit$rank)])
  x[, valid_columns, drop = FALSE]
}

empty_terms <- function() {
  data.frame(
    Group = character(), Group_Label = character(),
    Term_Internal = character(), Dietary_Variable = character(),
    F = numeric(), p_value = numeric(), Delta_R2 = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_scores <- function() {
  data.frame(
    Group = character(), Group_Label = character(),
    SEQN = character(), CAP1 = numeric(), CAP2 = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_biplot <- function() {
  data.frame(
    Group = character(), Group_Label = character(),
    Term_Internal = character(), Dietary_Variable = character(),
    CAP1 = numeric(), CAP2 = numeric(),
    stringsAsFactors = FALSE
  )
}

# Extract the first two ordination coordinates, allowing a single axis.
two_axes <- function(mod, display, row_ids) {
  s <- tryCatch(
    as.matrix(vegan::scores(mod, display = display,
                            choices = 1:2, scaling = 1)),
    error = function(e) NULL
  )
  if (is.null(s)) {
    s <- matrix(NA_real_, nrow = length(row_ids), ncol = 2,
                dimnames = list(row_ids, c("CAP1", "CAP2")))
  } else {
    if (ncol(s) < 2) s <- cbind(s, CAP2 = NA_real_)
    s <- s[, 1:2, drop = FALSE]
    colnames(s) <- c("CAP1", "CAP2")
  }
  s
}

run_dbrda_for_group <- function(group_code, seed) {
  set.seed(seed)
  group_label <- unname(group_labels[group_code])

  idx <- if (group_code == "Overall") {
    seq_len(nrow(df))
  } else {
    which(as.character(df[[group_var]]) == group_code)
  }
  n <- length(idx)

  blank_result <- function(note) {
    list(
      summary = data.frame(
        Group = group_code, Group_Label = group_label, n = n,
        R2 = NA_real_, R2_adj = NA_real_,
        F_overall = NA_real_, p_overall = NA_real_,
        constrained_axes = NA_integer_, Note = note
      ),
      sequential = empty_terms(), marginal = empty_terms(),
      sites = empty_scores(), biplot = empty_biplot()
    )
  }

  if (n < 4) return(blank_result("Insufficient samples"))

  comm_g <- clr_all[idx, , drop = FALSE]
  env_g <- diet_z[idx, , drop = FALSE]

  # Remove group-specific zero-variance predictors and exact dependencies.
  sd_g <- vapply(env_g, sd, numeric(1), na.rm = TRUE)
  env_g <- env_g[, is.finite(sd_g) & sd_g > 0, drop = FALSE]
  env_g <- drop_linear_dependencies(env_g)
  if (ncol(env_g) == 0) {
    return(blank_result("No valid dietary predictors"))
  }

  # All remaining dietary variables enter the same constrained model.
  mod <- vegan::capscale(
    comm_g ~ ., data = as.data.frame(env_g), distance = "euclidean"
  )

  model_r2 <- vegan::RsquareAdj(mod)
  r2_full <- as.numeric(model_r2$r.squared)
  r2_adj <- as.numeric(model_r2$adj.r.squared)

  # Test the overall constrained model.
  test_overall <- vegan::anova.cca(mod, permutations = n_perm)

  # Two variable-specific tests on the SAME fitted model:
  # (a) sequential: order-dependent; (b) marginal: adjusted for other terms.
  test_sequential <- vegan::anova.cca(
    mod, permutations = n_perm, by = "terms"
  )
  test_marginal <- vegan::anova.cca(
    mod, permutations = n_perm, by = "margin"
  )

  # Incremental (semi-partial) R2: full R2 minus leave-one-term-out R2.
  delta_r2 <- vapply(names(env_g), function(term) {
    remaining <- setdiff(names(env_g), term)
    if (length(remaining) == 0) return(r2_full)
    mod_reduced <- vegan::capscale(
      comm_g ~ .,
      data = as.data.frame(env_g[, remaining, drop = FALSE]),
      distance = "euclidean"
    )
    r2_reduced <- as.numeric(vegan::RsquareAdj(mod_reduced)$r.squared)
    if (!is.finite(r2_reduced)) return(NA_real_)
    max(0, r2_full - r2_reduced)
  }, numeric(1))

  format_test <- function(test) {
    tab <- as.data.frame(test)
    tab$Term_Internal <- rownames(tab)
    tab <- tab[tab$Term_Internal %in% names(env_g), , drop = FALSE]
    original <- diet_map$Dietary_Variable[
      match(tab$Term_Internal, diet_map$Internal)
    ]
    data.frame(
      Group = group_code,
      Group_Label = group_label,
      Term_Internal = tab$Term_Internal,
      Dietary_Variable = original,
      F = as.numeric(tab$F),
      p_value = as.numeric(tab[["Pr(>F)"]]),
      Delta_R2 = unname(delta_r2[tab$Term_Internal]),
      stringsAsFactors = FALSE
    )
  }

  # Scores are exported as numerical tables; no plotting is performed.
  site_sc <- two_axes(mod, "sites", rownames(comm_g))
  sites_tbl <- data.frame(
    Group = group_code, Group_Label = group_label,
    SEQN = rownames(site_sc),
    CAP1 = site_sc[, 1], CAP2 = site_sc[, 2],
    stringsAsFactors = FALSE
  )

  bp_sc <- two_axes(mod, "bp", names(env_g))
  biplot_tbl <- data.frame(
    Group = group_code, Group_Label = group_label,
    Term_Internal = rownames(bp_sc),
    Dietary_Variable = diet_map$Dietary_Variable[
      match(rownames(bp_sc), diet_map$Internal)
    ],
    CAP1 = bp_sc[, 1], CAP2 = bp_sc[, 2],
    stringsAsFactors = FALSE
  )

  list(
    summary = data.frame(
      Group = group_code, Group_Label = group_label, n = n,
      R2 = r2_full, R2_adj = r2_adj,
      F_overall = as.numeric(test_overall$F[1]),
      p_overall = as.numeric(test_overall[["Pr(>F)"]][1]),
      constrained_axes = length(vegan::eigenvals(mod, constrained = TRUE)),
      Note = NA_character_, stringsAsFactors = FALSE
    ),
    sequential = format_test(test_sequential),
    marginal = format_test(test_marginal),
    sites = sites_tbl,
    biplot = biplot_tbl
  )
}

# 4. Run the overall and oral-condition-specific models -------

core_count <- parallel::detectCores()
n_cores <- if (is.na(core_count)) 1L else {
  max(1L, min(length(group_codes), core_count - 1L))
}

run_one <- function(i) {
  run_dbrda_for_group(group_codes[i], seed = random_seed + i)
}

if (.Platform$OS.type == "windows" || n_cores == 1L) {
  results <- lapply(seq_along(group_codes), run_one)
} else {
  results <- parallel::mclapply(
    seq_along(group_codes), run_one,
    mc.cores = n_cores, mc.set.seed = FALSE
  )
}

summary_all <- do.call(rbind, lapply(results, `[[`, "summary"))
sequential_all <- do.call(rbind, lapply(results, `[[`, "sequential"))
marginal_all <- do.call(rbind, lapply(results, `[[`, "marginal"))

# 5. Export all statistical results to ONE Excel workbook -----

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "dbRDA Summary")
openxlsx::writeData(wb, "dbRDA Summary", summary_all)

openxlsx::addWorksheet(wb, "Terms Sequential")
openxlsx::writeData(wb, "Terms Sequential", sequential_all)

openxlsx::addWorksheet(wb, "Terms Marginal")
openxlsx::writeData(wb, "Terms Marginal", marginal_all)

for (i in seq_along(group_codes)) {
  label <- unname(group_labels[group_codes[i]])
  sites_sheet <- paste0("Sites_", label)
  biplot_sheet <- paste0("Biplot_", label)

  openxlsx::addWorksheet(wb, sites_sheet)
  openxlsx::writeData(wb, sites_sheet, results[[i]]$sites)

  openxlsx::addWorksheet(wb, biplot_sheet)
  openxlsx::writeData(wb, biplot_sheet, results[[i]]$biplot)
}

openxlsx::addWorksheet(wb, "Dietary Columns Used")
openxlsx::writeData(wb, "Dietary Columns Used", diet_map)

meta <- data.frame(
  Setting = c(
    "Distance", "Energy adjustment", "Dietary standardization",
    "Groups", "Permutations per test", "Random seed base",
    "Variable P values", "Variable explained variance",
    "Total retained sample size", "vegan version"
  ),
  Value = c(
    "Aitchison (Euclidean distance on CLR)",
    "Dietary intake per 1,000 kcal",
    "Global z-score before subgroup analyses",
    paste(group_codes, collapse = ", "),
    as.character(n_perm), as.character(random_seed),
    "Sequential (by=terms) and marginal (by=margin)",
    "Unadjusted full R2 minus leave-one-variable-out R2",
    as.character(nrow(df)), as.character(packageVersion("vegan"))
  ),
  stringsAsFactors = FALSE
)
openxlsx::addWorksheet(wb, "Meta")
openxlsx::writeData(wb, "Meta", meta)

openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
