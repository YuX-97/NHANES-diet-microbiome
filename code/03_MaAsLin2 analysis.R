# ============================================================
# MaAsLin2 — Energy-adjusted dietary analysis
#
# Model:
#   Taxon_CLR ~ dietary predictor, adjusted for age + sex + oral condition
#
# Settings:
#   - All dietary variables are energy-adjusted per 1000 kcal
#   - FDR = MaAsLin2 qval (BH across taxa for each predictor)
#
# Inputs:
#   - CLR microbiome table
#   - Metadata table
#
# Output:
#   - Excel file with MaAsLin2 results for all energy-adjusted dietary predictors
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(janitor)
  library(openxlsx)
  library(Maaslin2)
  library(tibble)
})

microbiome_path <- "CLR_3254 for MaAsLin2.tsv"
metadata_path <- "metadata.xlsx"
metadata_sheet <- "MaAsLin_3254"

output_xlsx <- "MaAsLin2_Diet_EnergyAdjusted.xlsx"

food_cols_idx <- 18:59
energy_col <- "Energy_intake_kcal"

covars_cont <- c("RIDAGEYR")
covars_cat <- c("Oral_CDC", "RIAGENDR")

micro_df <- read_tsv(microbiome_path, show_col_types = FALSE)
meta_df <- read_xlsx(metadata_path, sheet = metadata_sheet)

common_ids <- intersect(micro_df$SEQN, meta_df$SEQN)
micro_df <- micro_df %>% filter(SEQN %in% common_ids)
meta_df <- meta_df %>% filter(SEQN %in% common_ids)

taxa_cols <- setdiff(colnames(micro_df), "SEQN")
micro_df[taxa_cols] <- lapply(micro_df[taxa_cols], as.numeric)

feat_mat <- micro_df %>%
  arrange(SEQN) %>%
  column_to_rownames("SEQN") %>%
  as.matrix() %>%
  t()

orig_taxa <- rownames(feat_mat)
safe_taxa <- make.names(orig_taxa, unique = TRUE)

taxa_map <- tibble(
  TaxaSafe = safe_taxa,
  Taxa = orig_taxa
)

rownames(feat_mat) <- safe_taxa
features <- feat_mat
taxa_all <- orig_taxa

meta_df <- meta_df %>%
  arrange(SEQN) %>%
  column_to_rownames("SEQN")

food_cols_no_seqn <- food_cols_idx - 1
food_names <- colnames(meta_df)[food_cols_no_seqn]
food_names <- setdiff(food_names, energy_col)

make_pseudo <- function(x) {
  paste0("diet__", make_clean_names(x))
}

food_pseudo <- setNames(make_pseudo(food_names), food_names)

diet_map <- tibble(
  Predictor_Original = names(food_pseudo),
  Predictor_Pseudo = unname(food_pseudo)
)

build_meta_ea <- function() {
  meta_use <- meta_df
  
  meta_use[covars_cont] <- lapply(meta_use[covars_cont], as.numeric)
  meta_use[covars_cat] <- lapply(meta_use[covars_cat], as.factor)
  
  energy_k <- suppressWarnings(as.numeric(meta_use[[energy_col]]))
  denom <- energy_k / 1000
  denom[denom <= 0 | !is.finite(denom)] <- NA_real_
  
  food_block <- meta_use[, food_names, drop = FALSE]
  food_block <- as.data.frame(lapply(food_block, function(x) suppressWarnings(as.numeric(x))))
  food_block <- as.data.frame(sweep(food_block, 1, denom, "/"))
  
  colnames(food_block) <- unname(food_pseudo[colnames(food_block)])
  
  meta_use <- meta_use %>%
    select(-all_of(food_names)) %>%
    cbind(food_block)
  
  meta_use
}

run_energy_adjusted_model <- function() {
  meta_use <- build_meta_ea()
  all_res <- list()
  
  for (pred in diet_map$Predictor_Pseudo) {
    message("Running predictor: ", pred)
    
    needed <- c(pred, covars_cont, covars_cat)
    keep <- complete.cases(meta_use[, needed, drop = FALSE])
    ids <- rownames(meta_use)[keep]
    
    if (length(ids) < 50) next
    
    feat_sub <- features[, ids, drop = FALSE]
    meta_sub <- meta_use[ids, , drop = FALSE]
    
    out_dir <- file.path(tempdir(), paste0("maaslin2_ea_", pred))
    dir.create(out_dir, showWarnings = FALSE)
    
    args <- list(
      input_data = feat_sub,
      input_metadata = meta_sub,
      output = out_dir,
      analysis_method = "LM",
      normalization = "NONE",
      transform = "NONE",
      correction = "BH",
      fixed_effects = c(pred, covars_cont, covars_cat),
      random_effects = NULL,
      min_prevalence = 0,
      min_abundance = 0,
      standardize = FALSE
    )
    
    supported <- names(formals(Maaslin2))
    if ("plot_heatmap" %in% supported) args$plot_heatmap <- FALSE
    if ("plot_scatter" %in% supported) args$plot_scatter <- FALSE
    if ("plot_associations" %in% supported) args$plot_associations <- FALSE
    
    invisible(capture.output({
      suppressMessages({
        fit <- do.call(Maaslin2, args)
      })
    }, type = "output"))
    
    res_file <- file.path(out_dir, "all_results.tsv")
    if (!file.exists(res_file)) next
    
    res_tbl <- read_tsv(res_file, show_col_types = FALSE)
    
    main <- res_tbl %>%
      filter(metadata == pred) %>%
      transmute(
        TaxaSafe = feature,
        beta = coef,
        pval = pval,
        FDR = qval
      ) %>%
      left_join(taxa_map, by = "TaxaSafe") %>%
      transmute(
        Taxa = Taxa,
        beta = beta,
        pval = pval,
        FDR = FDR
      )
    
    main$Predictor_Pseudo <- pred
    main$n_kept <- length(ids)
    
    all_res[[pred]] <- main
  }
  
  pred_order <- diet_map$Predictor_Original
  
  final <- bind_rows(all_res) %>%
    left_join(diet_map, by = "Predictor_Pseudo") %>%
    mutate(
      Predictor_Original = factor(Predictor_Original, levels = pred_order),
      Taxa = factor(Taxa, levels = taxa_all)
    ) %>%
    arrange(Predictor_Original, Taxa) %>%
    select(
      Predictor_Original,
      Predictor_Pseudo,
      Taxa,
      beta,
      pval,
      FDR,
      n_kept
    )
  
  wb <- createWorkbook()
  addWorksheet(wb, "MaAsLin2_Results")
  writeData(wb, "MaAsLin2_Results", final)
  
  addWorksheet(wb, "Metadata")
  meta_out <- tibble(
    microbiome_file = microbiome_path,
    metadata_file = metadata_path,
    metadata_sheet = metadata_sheet,
    dietary_columns = paste0(min(food_cols_idx), "-", max(food_cols_idx)),
    energy_adjustment = paste0("All dietary predictors adjusted per 1000 kcal using ", energy_col),
    continuous_covariates = paste(covars_cont, collapse = ", "),
    categorical_covariates = paste(covars_cat, collapse = ", "),
    normalization = "NONE",
    transform = "NONE",
    standardize = FALSE,
    analysis_method = "LM"
  )
  writeData(wb, "Metadata", meta_out)
  
  addWorksheet(wb, "Predictor_Map")
  writeData(wb, "Predictor_Map", diet_map)
  
  addWorksheet(wb, "Taxa_Map")
  writeData(wb, "Taxa_Map", taxa_map)
  
  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
}

run_energy_adjusted_model()