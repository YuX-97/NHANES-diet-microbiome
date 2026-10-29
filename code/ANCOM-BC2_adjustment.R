# ============================================================
# ANCOM-BC2: genus-level differential abundance by oral condition
# Unadjusted and Model 1 (seven covariates)
# Reference: Oral_CDC = 1 (Health); contrasts: 2/3/4 vs 1
# ============================================================

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(phyloseq)
  library(openxlsx)
})

# 1. Input and analysis variables -----------------------------

input_file <- "genus_read_counts_with_covariates.tsv"
output_file <- "ANCOMBC2_unadjusted_and_model1.xlsx"

id_var <- "SEQN"
group_var <- "Oral_CDC"
reference_level <- "1"
taxa_vars <- character(0)  # Supply genus read-count variable NAMES.

continuous_covars <- c("RIDAGEYR", "Countall_teeth", "INDFMPIR")
categorical_covars <- c(
  "RIAGENDR", "SMD_status_1", "RIDRETH1_co", "DMDEDUC2_co"
)
model1_covars <- c(continuous_covars, categorical_covars)

# 2. Import counts and metadata -------------------------------

dat <- read.delim(
  input_file, header = TRUE, sep = "\t", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

stopifnot(length(taxa_vars) > 0)
required_vars <- c(id_var, group_var, taxa_vars, model1_covars)
stopifnot(all(required_vars %in% names(dat)))

# Restrict to the four study groups; use SEQN as the sample identifier.
dat <- dat[as.character(dat[[group_var]]) %in% c("1", "2", "3", "4"), , drop = FALSE]
stopifnot(nrow(dat) > 0, !anyNA(dat[[id_var]]), !anyDuplicated(dat[[id_var]]))
rownames(dat) <- as.character(dat[[id_var]])

count_data <- as.data.frame(lapply(
  dat[, taxa_vars, drop = FALSE],
  function(x) suppressWarnings(as.numeric(as.character(x)))
), check.names = FALSE)
names(count_data) <- taxa_vars
rownames(count_data) <- rownames(dat)

# Preserve the original script's treatment of missing read counts.
count_data[is.na(count_data)] <- 0
stopifnot(all(is.finite(as.matrix(count_data))), all(as.matrix(count_data) >= 0))

prepare_metadata <- function(adjusted = FALSE) {
  vars <- c(group_var, if (adjusted) model1_covars)
  meta <- dat[, vars, drop = FALSE]

  meta[[group_var]] <- factor(as.character(meta[[group_var]]),
                              levels = c("1", "2", "3", "4"))
  if (adjusted) {
    for (v in continuous_covars) {
      meta[[v]] <- suppressWarnings(as.numeric(as.character(meta[[v]])))
    }
    for (v in categorical_covars) {
      meta[[v]] <- factor(as.character(meta[[v]]))
    }
  }

  meta <- meta[complete.cases(meta), , drop = FALSE]
  meta[[group_var]] <- droplevels(meta[[group_var]])
  stopifnot(reference_level %in% levels(meta[[group_var]]))
  meta[[group_var]] <- relevel(meta[[group_var]], ref = reference_level)
  meta
}

# 3. Fit ANCOM-BC2 --------------------------------------------

find_result_column <- function(result, prefix, term) {
  candidates <- c(paste0(prefix, term), paste0(prefix, make.names(term)))
  hit <- candidates[candidates %in% names(result)]
  if (length(hit)) return(hit[1])
  hit <- grep(paste0("^", prefix, ".*", term, "$"),
              names(result), value = TRUE)
  if (length(hit)) return(hit[1])
  NA_character_
}

run_model <- function(adjusted = FALSE) {
  meta <- prepare_metadata(adjusted)
  counts <- t(as.matrix(count_data[rownames(meta), , drop = FALSE]))
  storage.mode(counts) <- "numeric"
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  stopifnot(nrow(counts) > 0)

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::sample_data(meta)
  )

  fix_formula <- if (adjusted) {
    paste(c(model1_covars, group_var), collapse = " + ")
  } else {
    group_var
  }

  set.seed(123)
  fit <- ANCOMBC::ancombc2(
    data = ps,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = fix_formula,
    rand_formula = NULL,
    p_adj_method = "BH",
    pseudo_sens = FALSE,
    prv_cut = 0,
    lib_cut = 0,
    s0_perc = 0.05,
    group = group_var,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 1,
    verbose = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    iter_control = list(tol = 1e-5, max_iter = 100, verbose = FALSE),
    em_control = list(tol = 1e-5, max_iter = 100),
    lme_control = NULL,
    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
    trend_control = NULL
  )

  result <- as.data.frame(fit$res)
  stopifnot("taxon" %in% names(result))

  tables <- list()
  for (level in c("2", "3", "4")) {
    term <- paste0(group_var, level)
    lfc_col <- find_result_column(result, "lfc_", term)
    p_col <- find_result_column(result, "p_", term)
    q_col <- find_result_column(result, "q_", term)
    if (is.na(p_col)) p_col <- find_result_column(result, "p_val_", term)
    if (is.na(q_col)) q_col <- find_result_column(result, "q_val_", term)
    if (anyNA(c(lfc_col, p_col, q_col))) {
      warning("ANCOM-BC2 coefficient columns not found for ", level, " vs 1")
      next
    }

    table_name <- paste0(if (adjusted) "Model1_" else "Unadjusted_",
                         level, "_vs_1")
    tables[[table_name]] <- data.frame(
      Genus = result$taxon,
      log2FoldChange = result[[lfc_col]] / log(2),
      pvalue = result[[p_col]],
      padj = result[[q_col]],
      stringsAsFactors = FALSE
    )
  }

  list(
    tables = tables,
    meta = data.frame(
      Model = if (adjusted) "Model 1" else "Unadjusted",
      Formula = fix_formula,
      n_samples = phyloseq::nsamples(ps),
      n_taxa = phyloseq::ntaxa(ps),
      stringsAsFactors = FALSE
    )
  )
}

unadjusted <- run_model(adjusted = FALSE)
model1 <- run_model(adjusted = TRUE)

# 4. Export both models to one Excel workbook -----------------

results <- c(unadjusted$tables, model1$tables)
results[["Model Information"]] <- rbind(unadjusted$meta, model1$meta)
openxlsx::write.xlsx(results, file = output_file, overwrite = TRUE)
