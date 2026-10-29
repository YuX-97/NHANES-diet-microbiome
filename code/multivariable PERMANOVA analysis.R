# ============================================================
# Multivariable PERMANOVA: marginal (adjusted) effects
# Aitchison distance = Euclidean distance on CLR-transformed taxa
# The dietary intake variables are entered as one matrix/block term.
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(openxlsx)
})

# 1. Input, output, and analysis settings ----------------------

input_file <- "analysis_data_with_CLR.tsv"
output_file <- "PERMANOVA_multivariable_marginal.xlsx"

n_perm <- 1001
random_seed <- 3000
available_cores <- parallel::detectCores()
n_cores <- if (is.na(available_cores)) 1L else max(1L, available_cores - 1L)

# Specify column NAMES (not positions) from the analysis dataset.
# These vectors are intentionally left as placeholders for sharing.
cat_vars <- character(0)       # Categorical covariate names
cont_vars <- character(0)      # Continuous covariate names
diet_vars <- character(0)      # Dietary intake variable names
clr_taxa_vars <- character(0)  # CLR-transformed genus variable names
energy_var <- "Energy_intake_kcal"

# 2. Import data and prepare variables ------------------------

dat <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

stopifnot(
  length(cat_vars) > 0,
  length(cont_vars) > 0,
  length(diet_vars) > 0,
  length(clr_taxa_vars) > 0
)

required_vars <- c(cat_vars, cont_vars, diet_vars, clr_taxa_vars, energy_var)
stopifnot(all(required_vars %in% names(dat)))

# CLR taxa matrix used to calculate Aitchison distance.
clr_all <- as.matrix(dat[, clr_taxa_vars, drop = FALSE])
storage.mode(clr_all) <- "numeric"

# Categorical and continuous covariates.
cat_df <- as.data.frame(lapply(dat[, cat_vars, drop = FALSE], as.factor))
cont_df <- as.data.frame(lapply(
  dat[, cont_vars, drop = FALSE],
  function(x) suppressWarnings(as.numeric(x))
))

cat_internal <- paste0("CAT_", seq_along(cat_vars))
cont_internal <- paste0("CONT_", seq_along(cont_vars))
names(cat_df) <- cat_internal
names(cont_df) <- cont_internal
other_df <- data.frame(cat_df, cont_df, check.names = FALSE)

term_map <- data.frame(
  term_internal = c(cat_internal, cont_internal, "Dietary_intake_block_EA"),
  variable = c(cat_vars, cont_vars, "Dietary intake block (EA, per 1000 kcal)"),
  var_type = c(
    rep("categorical", length(cat_vars)),
    rep("continuous", length(cont_vars)),
    "dietary_block"
  ),
  stringsAsFactors = FALSE
)

# 3. Energy-adjust dietary intake per 1000 kcal ---------------

energy <- suppressWarnings(as.numeric(dat[[energy_var]]))
energy[energy <= 0] <- NA_real_

diet_raw <- dat[, diet_vars, drop = FALSE]
diet_ea_all <- as.data.frame(lapply(diet_raw, function(x) {
  suppressWarnings(as.numeric(x)) / (energy / 1000)
}))

diet_internal <- paste0("DIET_", seq_along(diet_vars))
names(diet_ea_all) <- diet_internal
diet_map <- data.frame(
  diet_internal = diet_internal,
  diet_original = diet_vars,
  stringsAsFactors = FALSE
)

# 4. Complete-case filtering and variable checks --------------

keep <- complete.cases(clr_all) &
  complete.cases(other_df) &
  !is.na(energy) &
  complete.cases(diet_ea_all)

stopifnot(sum(keep) >= 20)

clr_use <- clr_all[keep, , drop = FALSE]
other_use <- other_df[keep, , drop = FALSE]
diet_use <- diet_ea_all[keep, , drop = FALSE]

for (v in cat_internal) {
  other_use[[v]] <- droplevels(other_use[[v]])
}

# Remove predictors that have no variation in the analysis sample.
cat_valid <- vapply(cat_internal, function(v) {
  nlevels(other_use[[v]]) >= 2
}, logical(1))
cont_valid <- vapply(cont_internal, function(v) {
  sd(other_use[[v]], na.rm = TRUE) > 0
}, logical(1))
diet_valid <- vapply(diet_use, function(x) {
  sd(x, na.rm = TRUE) > 0
}, logical(1))

valid_other_terms <- c(cat_internal[cat_valid], cont_internal[cont_valid])
diet_use <- diet_use[, diet_valid, drop = FALSE]
stopifnot(length(valid_other_terms) > 0, ncol(diet_use) > 0)

term_map_used <- term_map[
  term_map$term_internal %in% c(valid_other_terms, "Dietary_intake_block_EA"),
  , drop = FALSE
]
diet_map_used <- diet_map[
  diet_map$diet_internal %in% names(diet_use),
  , drop = FALSE
]

# 5. Build the multivariable model ----------------------------

model_df <- other_use[, valid_other_terms, drop = FALSE]
model_df$Dietary_intake_block_EA <- I(as.matrix(diet_use))

aitchison_distance <- dist(clr_use, method = "euclidean")
model_terms <- c(valid_other_terms, "Dietary_intake_block_EA")
model_formula <- as.formula(paste(
  "aitchison_distance ~", paste(model_terms, collapse = " + ")
))

# 6. Marginal PERMANOVA (adjusted for all other terms) --------

set.seed(random_seed)
fit_margin <- vegan::adonis2(
  model_formula,
  data = model_df,
  permutations = n_perm,
  by = "margin",
  parallel = n_cores
)

# 7. Format and export marginal results -----------------------

fit_df <- as.data.frame(fit_margin)
fit_df$term_internal <- rownames(fit_df)
fit_df <- fit_df[
  !fit_df$term_internal %in% c("Residual", "Total", "Model"),
  , drop = FALSE
]

term_idx <- match(fit_df$term_internal, term_map_used$term_internal)
stopifnot(!anyNA(term_idx))

res_margin <- data.frame(
  variable = term_map_used$variable[term_idx],
  term_internal = fit_df$term_internal,
  var_type = term_map_used$var_type[term_idx],
  n_used = nrow(model_df),
  df = as.numeric(fit_df$Df),
  R2 = as.numeric(fit_df$R2),
  F = as.numeric(fit_df$F),
  p_value = as.numeric(fit_df[["Pr(>F)"]]),
  stringsAsFactors = FALSE
)
res_margin <- res_margin[order(res_margin$R2, decreasing = TRUE), ]

meta <- data.frame(
  setting = c(
    "distance", "test", "permutations", "random_seed",
    "n_samples_used", "n_diet_variables_used", "energy_adjustment"
  ),
  value = c(
    "Aitchison (Euclidean on CLR)", "adonis2: by = margin",
    as.character(n_perm), as.character(random_seed),
    as.character(nrow(model_df)), as.character(ncol(diet_use)),
    "Dietary variables per 1000 kcal; one matrix/block term"
  ),
  stringsAsFactors = FALSE
)

openxlsx::write.xlsx(
  list(
    Multivariable_margin_R2 = res_margin,
    Variables_used = term_map_used,
    Diet_variables_used = diet_map_used,
    Meta = meta
  ),
  file = output_file,
  overwrite = TRUE
)
