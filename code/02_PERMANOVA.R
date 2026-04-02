# ============================================================
# NHANES oral rinse microbiome — Unadjusted PERMANOVA R2 scan
# Parallel version (macOS/Linux via mclapply)
#
# Input: CLR-transformed taxa table included in the same TSV
# Distance: Aitchison = Euclidean distance on CLR
#
# Analyses:
# 1) Single-variable PERMANOVA (unadjusted/total R2)
#    - categorical predictors: columns 1-9
#    - continuous predictors : columns 11-15
# 2) Dietary intake block PERMANOVA
#    - diet variables: columns 16-53
#    - energy adjustment: per 1000 kcal using Energy_intake_kcal
#    - model-level (block) R2 and p
#
# Output: Excel with results + metadata (+ optional Diet_EA_used sheet)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(vegan)
  library(openxlsx)
  library(parallel)
})

input_path <- "CLR_3254 for PERMANOVA.tsv"
output_xlsx <- sub("\\.tsv$", "_PERMANOVA_Unadjusted_R2.xlsx", input_path)

if (identical(output_xlsx, input_path)) {
  output_xlsx <- "PERMANOVA_Unadjusted_R2.xlsx"
}

n_perm <- 1001
n_cores <- max(1, parallel::detectCores() - 1)

cat_cols_idx <- 1:9
cont_cols_idx <- 11:15
energy_col_name <- "Energy_intake_kcal"
diet_cols_idx <- 16:53
taxa_cols_idx <- 54:130

export_diet_ea_sheet <- TRUE

dat <- read.delim(
  input_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

stopifnot(ncol(dat) >= max(taxa_cols_idx))
stopifnot(energy_col_name %in% colnames(dat))

clr_mat_all <- dat[, taxa_cols_idx, drop = FALSE]
clr_mat_all <- apply(clr_mat_all, 2, function(x) suppressWarnings(as.numeric(x)))
clr_mat_all <- as.matrix(clr_mat_all)

run_single_permanova <- function(dat_in, clr_mat_in, var_name,
                                 var_type = c("categorical", "continuous"),
                                 n_perm = 1001,
                                 seed = 1) {
  var_type <- match.arg(var_type)
  
  if (!var_name %in% colnames(dat_in)) {
    return(tibble(
      variable = var_name,
      var_type = var_type,
      n_used = NA_integer_,
      df = NA_real_,
      R2 = NA_real_,
      F = NA_real_,
      p_value = NA_real_,
      note = "variable_not_found"
    ))
  }
  
  x <- dat_in[[var_name]]
  
  if (var_type == "categorical") {
    x <- as.factor(x)
  } else {
    x <- suppressWarnings(as.numeric(x))
  }
  
  keep <- !is.na(x) & complete.cases(clr_mat_in)
  
  if (sum(keep) < 5) {
    return(tibble(
      variable = var_name,
      var_type = var_type,
      n_used = sum(keep),
      df = NA_real_,
      R2 = NA_real_,
      F = NA_real_,
      p_value = NA_real_,
      note = "too_few_samples_after_filtering"
    ))
  }
  
  if (var_type == "categorical") {
    x2 <- droplevels(x[keep])
    if (nlevels(x2) < 2) {
      return(tibble(
        variable = var_name,
        var_type = var_type,
        n_used = sum(keep),
        df = NA_real_,
        R2 = NA_real_,
        F = NA_real_,
        p_value = NA_real_,
        note = "fewer_than_2_levels_after_filtering"
      ))
    }
    x_use <- x2
  } else {
    x_use <- x[keep]
    if (sd(x_use, na.rm = TRUE) == 0) {
      return(tibble(
        variable = var_name,
        var_type = var_type,
        n_used = sum(keep),
        df = NA_real_,
        R2 = NA_real_,
        F = NA_real_,
        p_value = NA_real_,
        note = "zero_variance_after_filtering"
      ))
    }
  }
  
  clr_use <- clr_mat_in[keep, , drop = FALSE]
  d <- dist(clr_use, method = "euclidean")
  df_model <- data.frame(x = x_use)
  
  set.seed(seed)
  fit <- adonis2(d ~ x, data = df_model, permutations = n_perm, by = "terms")
  term_row <- fit[rownames(fit) == "x", , drop = FALSE]
  
  tibble(
    variable = var_name,
    var_type = var_type,
    n_used = nrow(df_model),
    df = as.numeric(term_row$Df),
    R2 = as.numeric(term_row$R2),
    F = as.numeric(term_row$F),
    p_value = as.numeric(term_row$`Pr(>F)`),
    note = ""
  )
}

run_diet_block <- function(dat_in, clr_mat_in,
                           energy_col = "Energy_intake_kcal",
                           diet_cols_idx = 16:53,
                           n_perm = 1001,
                           seed = 3000) {
  diet_names <- colnames(dat_in)[diet_cols_idx]
  
  if (energy_col %in% diet_names) {
    diet_names <- setdiff(diet_names, energy_col)
  }
  
  e <- suppressWarnings(as.numeric(dat_in[[energy_col]]))
  e[e <= 0] <- NA_real_
  
  diet_raw <- dat_in[, diet_names, drop = FALSE]
  diet_raw_num <- as.data.frame(lapply(diet_raw, function(v) suppressWarnings(as.numeric(v))))
  
  e_scale <- e / 1000
  diet_ea <- as.data.frame(lapply(diet_raw_num, function(v) v / e_scale))
  
  keep <- complete.cases(clr_mat_in) & !is.na(e) & complete.cases(diet_ea)
  
  if (sum(keep) < 20) {
    return(list(
      overall = tibble(
        variable = "Dietary_Intake_Block_EnergyAdjusted",
        var_type = "block_multivariable",
        n_used = sum(keep),
        df = NA_real_,
        R2 = NA_real_,
        F = NA_real_,
        p_value = NA_real_,
        note = "too_few_samples_after_filtering"
      ),
      diet_ea = NULL
    ))
  }
  
  diet_use <- diet_ea[keep, , drop = FALSE]
  nzv <- vapply(diet_use, function(v) sd(v, na.rm = TRUE) > 0, logical(1))
  diet_use <- diet_use[, nzv, drop = FALSE]
  
  if (ncol(diet_use) < 1) {
    return(list(
      overall = tibble(
        variable = "Dietary_Intake_Block_EnergyAdjusted",
        var_type = "block_multivariable",
        n_used = sum(keep),
        df = NA_real_,
        R2 = NA_real_,
        F = NA_real_,
        p_value = NA_real_,
        note = "all_diet_vars_zero_variance_after_filtering"
      ),
      diet_ea = NULL
    ))
  }
  
  clr_use <- clr_mat_in[keep, , drop = FALSE]
  d <- dist(clr_use, method = "euclidean")
  
  set.seed(seed)
  fit_model <- adonis2(d ~ ., data = diet_use, permutations = n_perm, by = NULL)
  model_row <- fit_model[rownames(fit_model) == "Model", , drop = FALSE]
  
  overall <- tibble(
    variable = "Dietary_Intake_Block_EnergyAdjusted",
    var_type = "block_multivariable",
    n_used = sum(keep),
    df = as.numeric(model_row$Df),
    R2 = as.numeric(model_row$R2),
    F = as.numeric(model_row$F),
    p_value = as.numeric(model_row$`Pr(>F)`),
    note = paste0(
      "Diet vars used: ", ncol(diet_use),
      " (from ", length(diet_names), "); energy-adjusted per 1000 kcal using ", energy_col
    )
  )
  
  list(overall = overall, diet_ea = diet_use)
}

cat_vars <- colnames(dat)[cat_cols_idx]
cont_vars <- colnames(dat)[cont_cols_idx]

cat_seeds <- 1000 + seq_along(cat_vars)
cont_seeds <- 2000 + seq_along(cont_vars)

res_cat <- bind_rows(parallel::mclapply(seq_along(cat_vars), function(i) {
  v <- cat_vars[i]
  run_single_permanova(
    dat,
    clr_mat_all,
    v,
    var_type = "categorical",
    n_perm = n_perm,
    seed = cat_seeds[i]
  )
}, mc.cores = n_cores))

res_cont <- bind_rows(parallel::mclapply(seq_along(cont_vars), function(i) {
  v <- cont_vars[i]
  run_single_permanova(
    dat,
    clr_mat_all,
    v,
    var_type = "continuous",
    n_perm = n_perm,
    seed = cont_seeds[i]
  )
}, mc.cores = n_cores))

diet_block <- run_diet_block(
  dat,
  clr_mat_all,
  energy_col = energy_col_name,
  diet_cols_idx = diet_cols_idx,
  n_perm = n_perm,
  seed = 3000
)

res_block <- diet_block$overall

res_all <- bind_rows(res_cat, res_cont, res_block) %>%
  mutate(
    p_value = as.numeric(p_value),
    R2 = as.numeric(R2)
  ) %>%
  arrange(desc(R2))

wb <- createWorkbook()

addWorksheet(wb, "PERMANOVA_R2")
writeData(wb, "PERMANOVA_R2", res_all)

addWorksheet(wb, "Metadata")
meta <- tibble(
  input_file = input_path,
  distance = "Aitchison_Euclidean_on_CLR",
  taxa_columns = paste0(min(taxa_cols_idx), "-", max(taxa_cols_idx)),
  categorical_predictors = paste(cat_vars, collapse = ", "),
  continuous_predictors = paste(cont_vars, collapse = ", "),
  diet_block = paste0(
    "Columns ", min(diet_cols_idx), "-", max(diet_cols_idx),
    "; energy-adjusted per 1000 kcal using ", energy_col_name
  ),
  permutations = n_perm,
  parallel_cores = n_cores
)
writeData(wb, "Metadata", meta)

if (export_diet_ea_sheet && !is.null(diet_block$diet_ea)) {
  addWorksheet(wb, "Diet_EA_Used")
  writeData(wb, "Diet_EA_Used", diet_block$diet_ea)
}

saveWorkbook(wb, output_xlsx, overwrite = TRUE)