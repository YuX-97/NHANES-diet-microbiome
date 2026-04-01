# ============================================================
# Multiple Parallel Mediation Analysis
#
# Predictor:
#   Added_sugar
#
# Outcome:
#   Oral_CDC
#   1 = oral health
#   2 = caries-only
#   3 = periodontitis
#   4 = caries+periodontitis
#
# Outcome contrasts analyzed:
#   2 vs 1
#   3 vs 1
#   4 vs 1
#
# Mediators:
#   cluster_1_score, cluster_2_score, cluster_3_score, cluster_4_score
#   from microbiome_cluster_scores.xlsx, sheet = "MaAsLin2_Z_mean"
#
# Estimation:
#   lavaan SEM with bootstrap confidence intervals for indirect,
#   direct, and total effects
#
# Parallel bootstrap:
#   multicore on macOS/Linux, snow on Windows
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(lavaan)
  library(parallel)
})

metadata_path <- "metadata.xlsx"
metadata_sheet <- "Mediation"

med_path <- "microbiome_cluster_scores.xlsx"
output_xlsx <- "Mediation_AddedSugar_Microbiome_OralCDC.xlsx"

mediator_sheet_name <- "MaAsLin2_Z_mean"

all_cores <- parallel::detectCores(logical = TRUE)
n_cores <- max(1, all_cores - 1)
lavaan_parallel <- if (.Platform$OS.type == "windows") "snow" else "multicore"

cat("Detected cores:", all_cores, "| Using", n_cores, "cores for lavaan bootstrap (", lavaan_parallel, ")\n")

metadata_dat <- readxl::read_excel(metadata_path, sheet = metadata_sheet)

required_cols_metadata <- c("SEQN", "Oral_CDC", "Added_sugar")
missing_metadata <- setdiff(required_cols_metadata, colnames(metadata_dat))
if (length(missing_metadata) > 0) {
  stop("Metadata file is missing required columns: ", paste(missing_metadata, collapse = ", "))
}

if (anyDuplicated(metadata_dat$SEQN) > 0) {
  warning("Duplicated SEQN values detected in metadata.")
}

all_med_sheets <- openxlsx::getSheetNames(med_path)
if (!(mediator_sheet_name %in% all_med_sheets)) {
  stop(
    "Mediator sheet not found: ", mediator_sheet_name,
    ". Available sheets: ", paste(all_med_sheets, collapse = ", ")
  )
}

med_dat_full <- openxlsx::read.xlsx(med_path, sheet = mediator_sheet_name)

required_cols_med <- c(
  "SEQN",
  "cluster_1_score",
  "cluster_2_score",
  "cluster_3_score",
  "cluster_4_score"
)
missing_med <- setdiff(required_cols_med, colnames(med_dat_full))
if (length(missing_med) > 0) {
  stop("Mediator data is missing required columns: ", paste(missing_med, collapse = ", "))
}

lav_model_template <- "
  cluster_1_score ~ a1*Added_sugar
  cluster_2_score ~ a2*Added_sugar
  cluster_3_score ~ a3*Added_sugar
  cluster_4_score ~ a4*Added_sugar

  Y ~ cprime*Added_sugar + b1*cluster_1_score + b2*cluster_2_score + b3*cluster_3_score + b4*cluster_4_score

  ind1 := a1*b1
  ind2 := a2*b2
  ind3 := a3*b3
  ind4 := a4*b4

  ind_total := ind1 + ind2 + ind3 + ind4
  total := cprime + ind_total
"

fit_parallel_mediation_contrast <- function(
    target_label,
    ref_label = 1,
    metadata_dat,
    med_dat,
    lav_model_template,
    boot_n = 2000
) {
  cat("\n============================================\n")
  cat("Processing Oral_CDC contrast:", target_label, "vs", ref_label, "\n")
  cat("============================================\n")
  
  dat_merged <- metadata_dat %>%
    dplyr::inner_join(med_dat, by = "SEQN")
  
  keep_labels <- c(ref_label, target_label)
  dat_bin <- dat_merged %>%
    dplyr::filter(Oral_CDC %in% keep_labels) %>%
    dplyr::mutate(Y = ifelse(Oral_CDC == target_label, 1, 0))
  
  dat_model <- dat_bin %>%
    dplyr::select(
      SEQN,
      Y,
      Added_sugar,
      cluster_1_score,
      cluster_2_score,
      cluster_3_score,
      cluster_4_score
    ) %>%
    dplyr::filter(!is.na(Added_sugar))
  
  dat_model <- stats::na.omit(dat_model)
  
  cat("N used =", nrow(dat_model), "\n")
  if (nrow(dat_model) < 50) {
    warning(
      "Sample size for Oral_CDC ", target_label, " vs ", ref_label,
      " is < 50. Results may be unstable."
    )
  }
  
  fit <- lavaan::sem(
    model = lav_model_template,
    data = dat_model,
    se = "bootstrap",
    bootstrap = boot_n,
    meanstructure = TRUE,
    parallel = lavaan_parallel,
    ncpus = n_cores
  )
  
  pe <- lavaan::parameterEstimates(
    fit,
    standardized = FALSE,
    ci = TRUE,
    level = 0.95,
    boot.ci.type = "perc"
  )
  
  target_labels <- c(
    paste0("a", 1:4),
    paste0("b", 1:4),
    "cprime",
    paste0("ind", 1:4),
    "ind_total",
    "total"
  )
  
  res <- pe %>%
    dplyr::filter(label %in% target_labels) %>%
    dplyr::mutate(
      effect_type = dplyr::case_when(
        label %in% paste0("a", 1:4) ~ "a (X -> Mi)",
        label %in% paste0("b", 1:4) ~ "b (Mi -> Y)",
        label == "cprime" ~ "direct (c_prime)",
        label %in% paste0("ind", 1:4) ~ "indirect (Mi)",
        label == "ind_total" ~ "total indirect",
        label == "total" ~ "total (c)",
        TRUE ~ ""
      ),
      mediator = dplyr::case_when(
        label %in% c("a1", "b1", "ind1") ~ "cluster_1_score",
        label %in% c("a2", "b2", "ind2") ~ "cluster_2_score",
        label %in% c("a3", "b3", "ind3") ~ "cluster_3_score",
        label %in% c("a4", "b4", "ind4") ~ "cluster_4_score",
        label %in% c("cprime", "ind_total", "total") ~ NA_character_,
        TRUE ~ NA_character_
      )
    )
  
  total_row <- res %>% dplyr::filter(label == "total")
  total_est <- if (nrow(total_row) > 0) total_row$est[1] else NA_real_
  eps <- 1e-8
  
  contrast_text <- paste0(
    "Oral_CDC: ", target_label, " vs ", ref_label,
    " (Y = 1 if ", target_label, ", 0 if ", ref_label, ")"
  )
  
  res <- res %>%
    dplyr::mutate(
      proportion_of_total = dplyr::case_when(
        effect_type %in% c("indirect (Mi)", "total indirect") &
          !is.na(total_est) & abs(total_est) > eps ~ est / total_est,
        TRUE ~ NA_real_
      ),
      outcome_contrast = contrast_text
    )
  
  res <- res %>%
    dplyr::mutate(
      effect_order = dplyr::case_when(
        effect_type == "a (X -> Mi)" ~ 1,
        effect_type == "b (Mi -> Y)" ~ 2,
        effect_type == "indirect (Mi)" ~ 3,
        effect_type == "total indirect" ~ 4,
        effect_type == "direct (c_prime)" ~ 5,
        effect_type == "total (c)" ~ 6,
        TRUE ~ 99
      ),
      mediator_order = dplyr::case_when(
        mediator == "cluster_1_score" ~ 1,
        mediator == "cluster_2_score" ~ 2,
        mediator == "cluster_3_score" ~ 3,
        mediator == "cluster_4_score" ~ 4,
        TRUE ~ 99
      )
    ) %>%
    dplyr::arrange(effect_order, mediator_order) %>%
    dplyr::select(-effect_order, -mediator_order) %>%
    dplyr::select(
      outcome_contrast,
      label,
      effect_type,
      mediator,
      lhs,
      op,
      rhs,
      est,
      se,
      z,
      pvalue,
      ci.lower,
      ci.upper,
      proportion_of_total
    )
  
  attr(res, "N_used") <- nrow(dat_model)
  
  res
}

contrast_list <- list(
  list(target = 2, ref = 1, sheet_name = "Caries_vs_OralHealth"),
  list(target = 3, ref = 1, sheet_name = "Periodontitis_vs_OralHealth"),
  list(target = 4, ref = 1, sheet_name = "CariesPeriodontitis_vs_OralHealth")
)

wb <- openxlsx::createWorkbook()

for (cc in contrast_list) {
  target_label <- cc$target
  ref_label <- cc$ref
  sheet_name <- cc$sheet_name
  
  res_contrast <- fit_parallel_mediation_contrast(
    target_label = target_label,
    ref_label = ref_label,
    metadata_dat = metadata_dat,
    med_dat = med_dat_full,
    lav_model_template = lav_model_template,
    boot_n = 2000
  )
  
  N_used <- attr(res_contrast, "N_used")
  cat("Contrast", target_label, "vs", ref_label, ": N_used =", N_used, "\n")
  
  res_contrast_out <- res_contrast %>%
    dplyr::mutate(N_used = N_used)
  
  openxlsx::addWorksheet(wb, sheetName = sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = res_contrast_out)
}

openxlsx::addWorksheet(wb, "Metadata")
metadata_out <- dplyr::tibble(
  metadata_file = metadata_path,
  metadata_sheet = metadata_sheet,
  mediator_file = med_path,
  mediator_sheet = mediator_sheet_name,
  predictor = "Added_sugar",
  outcome = "Oral_CDC",
  oral_cdc_mapping = "1=oral health; 2=caries-only; 3=periodontitis; 4=caries+periodontitis",
  contrasts = "2 vs 1; 3 vs 1; 4 vs 1",
  bootstrap = 2000,
  parallel_mode = lavaan_parallel,
  n_cores = n_cores
)
openxlsx::writeData(wb, sheet = "Metadata", x = metadata_out)

openxlsx::saveWorkbook(wb, file = output_xlsx, overwrite = TRUE)

cat(
  "\nAll multiple parallel mediation models (Added_sugar -> microbiome clusters -> Oral_CDC contrasts) have finished. Results written to:\n",
  output_xlsx,
  "\n"
)