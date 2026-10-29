# ============================================================
# MaAsLin2: dietary intake associations with genus-level CLR data
# EA scenario only: food intakes per 1,000 kcal; food count unchanged
# Model 1: age, sex, smoking, tooth number, race, education, income
# Model 2: Model 1 + oral condition
# Each model is exported to a separate Excel workbook.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(janitor)
  library(openxlsx)
  library(Maaslin2)
})

# 1. Inputs and variable names --------------------------------

microbiome_file <- "genus_CLR_with_SEQN.tsv"
metadata_file <- "dietary_intake_with_covariates.xlsx"
metadata_sheet <- "analysis_data"

output_model1 <- "MaAsLin2_Model1_EA.xlsx"
output_model2 <- "MaAsLin2_Model2_EA.xlsx"

id_var <- "SEQN"
group_var <- "Oral_CDC"
energy_var <- "Energy_intake_kcal"

# Supply food-intake variable NAMES (rather than column positions).
# Replace food_count_var with its exact column name if it differs.
food_intake_vars <- character(0)
food_count_var <- "number of foods reported"
food_vars <- c(food_intake_vars, food_count_var)

continuous_covars <- c("RIDAGEYR", "Countall_teeth", "INDFMPIR")
categorical_covars <- c(
  "RIAGENDR", "SMD_status_1", "RIDRETH1_co", "DMDEDUC2_co"
)
model1_covars <- c(continuous_covars, categorical_covars)
model2_covars <- c(model1_covars, group_var)

min_samples <- 50

# 2. Read data and match samples by SEQN ----------------------

micro_df <- as.data.frame(readr::read_tsv(
  microbiome_file, show_col_types = FALSE
))
meta_df <- as.data.frame(readxl::read_excel(
  metadata_file, sheet = metadata_sheet
))

stopifnot(length(food_intake_vars) > 0)
stopifnot(id_var %in% names(micro_df), id_var %in% names(meta_df))
required_meta <- c(food_vars, energy_var, model2_covars)
stopifnot(all(required_meta %in% names(meta_df)))

micro_df[[id_var]] <- as.character(micro_df[[id_var]])
meta_df[[id_var]] <- as.character(meta_df[[id_var]])
stopifnot(!anyDuplicated(micro_df[[id_var]]),
          !anyDuplicated(meta_df[[id_var]]))

# Restrict to the four prespecified oral-condition groups.
meta_df <- meta_df[
  as.character(meta_df[[group_var]]) %in% c("1", "2", "3", "4"),
  , drop = FALSE
]
common_ids <- intersect(micro_df[[id_var]], meta_df[[id_var]])
stopifnot(length(common_ids) > 0)

micro_df <- micro_df[match(common_ids, micro_df[[id_var]]), , drop = FALSE]
meta_df <- meta_df[match(common_ids, meta_df[[id_var]]), , drop = FALSE]
rownames(meta_df) <- common_ids

# All non-ID columns in the CLR input are genus-level features.
taxa_vars <- setdiff(names(micro_df), id_var)
stopifnot(length(taxa_vars) > 0)
clr_df <- as.data.frame(lapply(
  micro_df[, taxa_vars, drop = FALSE],
  function(x) suppressWarnings(as.numeric(x))
), check.names = FALSE)
names(clr_df) <- taxa_vars
rownames(clr_df) <- common_ids

# MaAsLin2 expects features as rows and samples as columns.
feature_mat <- t(as.matrix(clr_df))
taxa_map <- data.frame(
  TaxaSafe = make.names(rownames(feature_mat), unique = TRUE),
  Taxa = rownames(feature_mat),
  stringsAsFactors = FALSE
)
rownames(feature_mat) <- taxa_map$TaxaSafe

# 3. Energy adjustment and dietary z-score --------------------

meta_use <- meta_df
for (v in continuous_covars) {
  meta_use[[v]] <- suppressWarnings(as.numeric(as.character(meta_use[[v]])))
}
for (v in c(categorical_covars, group_var)) {
  meta_use[[v]] <- factor(as.character(meta_use[[v]]))
}

food_df <- as.data.frame(lapply(
  meta_use[, food_vars, drop = FALSE],
  function(x) suppressWarnings(as.numeric(as.character(x)))
), check.names = FALSE)
names(food_df) <- food_vars

energy <- suppressWarnings(as.numeric(as.character(meta_use[[energy_var]])))
energy[!is.finite(energy) | energy <= 0] <- NA_real_

# Energy adjustment applies to food intakes, NOT to food-count data.
for (v in food_intake_vars) {
  food_df[[v]] <- food_df[[v]] / (energy / 1000)
}

safe_zscore <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mu) / s
}

# Z-score every dietary predictor across all matched study participants.
# CLR features are not z-standardized a second time.
food_df[] <- lapply(food_df, safe_zscore)

food_pseudo <- setNames(
  paste0("diet__", janitor::make_clean_names(food_vars)), food_vars
)
diet_map <- data.frame(
  Predictor_Original = names(food_pseudo),
  Predictor_Pseudo = unname(food_pseudo),
  stringsAsFactors = FALSE
)
names(food_df) <- unname(food_pseudo[names(food_df)])

meta_use <- cbind(
  meta_use[, setdiff(names(meta_use), food_vars), drop = FALSE],
  food_df
)

# 4. Run one MaAsLin2 model for each dietary predictor ---------

run_model <- function(model_name, fixed_covars, output_file) {
  result_list <- list()

  for (pred in diet_map$Predictor_Pseudo) {
    # Predictor-specific complete-case analysis.
    needed <- c(pred, fixed_covars)
    keep <- complete.cases(meta_use[, needed, drop = FALSE])
    ids <- rownames(meta_use)[keep]
    if (length(ids) < min_samples) next

    feature_sub <- feature_mat[, ids, drop = FALSE]
    metadata_sub <- meta_use[ids, , drop = FALSE]

    # MaAsLin2 writes temporary intermediate results to its output folder.
    work_dir <- tempfile(pattern = paste0("maaslin2_", model_name, "_"))
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

    args <- list(
      input_data = feature_sub,
      input_metadata = metadata_sub,
      output = work_dir,
      analysis_method = "LM",
      normalization = "NONE",
      transform = "NONE",
      correction = "BH",
      fixed_effects = c(pred, fixed_covars),
      random_effects = NULL,
      min_prevalence = 0,
      min_abundance = 0,
      standardize = FALSE
    )

    supported <- names(formals(Maaslin2::Maaslin2))
    if ("plot_heatmap" %in% supported) args$plot_heatmap <- FALSE
    if ("plot_scatter" %in% supported) args$plot_scatter <- FALSE
    if ("plot_associations" %in% supported) args$plot_associations <- FALSE

    invisible(capture.output(
      suppressMessages(do.call(Maaslin2::Maaslin2, args)),
      type = "output"
    ))

    results_path <- file.path(work_dir, "all_results.tsv")
    if (!file.exists(results_path)) next
    result <- as.data.frame(readr::read_tsv(
      results_path, show_col_types = FALSE
    ))

    # Extract only the main effect of the current dietary variable.
    main <- result[result$metadata == pred, , drop = FALSE]
    if (!nrow(main)) next

    result_list[[pred]] <- data.frame(
      Predictor_Pseudo = pred,
      Taxa = taxa_map$Taxa[match(main$feature, taxa_map$TaxaSafe)],
      beta = main$coef,
      pval = main$pval,
      FDR = main$qval,
      n_kept = length(ids),
      stringsAsFactors = FALSE
    )
  }

  if (length(result_list)) {
    final <- dplyr::bind_rows(result_list)
    final$Predictor_Original <- diet_map$Predictor_Original[
      match(final$Predictor_Pseudo, diet_map$Predictor_Pseudo)
    ]
    final <- final[order(
      match(final$Predictor_Original, diet_map$Predictor_Original),
      match(final$Taxa, taxa_map$Taxa)
    ), c("Predictor_Original", "Predictor_Pseudo", "Taxa",
         "beta", "pval", "FDR", "n_kept")]
  } else {
    final <- data.frame(
      Predictor_Original = character(), Predictor_Pseudo = character(),
      Taxa = character(), beta = numeric(), pval = numeric(),
      FDR = numeric(), n_kept = integer()
    )
  }

  settings <- data.frame(
    Setting = c(
      "Model", "Fixed covariates", "Dietary scenario",
      "Food intake adjustment", "Number of foods reported",
      "Dietary z-score", "Microbial CLR z-score",
      "MaAsLin2 normalization", "MaAsLin2 transform",
      "Analysis", "Correction", "Minimum complete-case sample size"
    ),
    Value = c(
      model_name, paste(fixed_covars, collapse = ", "), "EA only",
      "Per 1,000 kcal", "No energy adjustment; z-score only",
      "TRUE", "FALSE", "NONE", "NONE", "LM", "BH",
      as.character(min_samples)
    ),
    stringsAsFactors = FALSE
  )

  openxlsx::write.xlsx(
    list(MaAsLin2_results = final,
         Run_settings = settings,
         Predictor_mapping = diet_map),
    file = output_file,
    overwrite = TRUE
  )
}

# 5. Fit both models and save separate Excel files ------------

run_model("Model1", model1_covars, output_model1)
run_model("Model2", model2_covars, output_model2)
