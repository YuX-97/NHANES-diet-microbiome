# ============================================================
# Dietary intake differences by oral condition: linear regression
# Unadjusted and Model 1 (seven covariates)
# Contrasts: Caries, Periodontitis, Coexistence vs Health
# Food intakes: per 1,000 kcal, then z-score; food count: z-score only
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
})

# 1. Input, output, and variable names -------------------------

input_file <- "dietary_intake_with_covariates.xlsx"
input_sheet <- "analysis_data"
output_file <- "Dietary_intake_oral_conditions_unadjusted_model1.xlsx"

group_var <- "Oral_CDC"
energy_var <- "Energy_intake_kcal"
reference_level <- "1"
target_levels <- c("2", "3", "4")

# Supply food-intake column NAMES from the analysis dataset.
# The food-count variable is not energy-adjusted; use its exact column name.
food_intake_vars <- character(0)
food_count_var <- "number of foods reported"
food_vars <- c(food_intake_vars, food_count_var)

continuous_covars <- c("Countall_teeth", "RIDAGEYR", "INDFMPIR")
categorical_covars <- c(
  "RIAGENDR", "SMD_status_1", "RIDRETH1_co", "DMDEDUC2_co"
)
model1_covars <- c(continuous_covars, categorical_covars)

# 2. Read data and prepare dietary outcomes -------------------

dat <- as.data.frame(readxl::read_excel(input_file, sheet = input_sheet))
stopifnot(length(food_intake_vars) > 0)
required_vars <- c(group_var, energy_var, food_vars, model1_covars)
stopifnot(all(required_vars %in% names(dat)))

dat <- dat[as.character(dat[[group_var]]) %in% c("1", "2", "3", "4"),
           , drop = FALSE]

food_df <- as.data.frame(lapply(
  dat[, food_vars, drop = FALSE],
  function(x) suppressWarnings(as.numeric(as.character(x)))
), check.names = FALSE)
names(food_df) <- food_vars

energy <- suppressWarnings(as.numeric(as.character(dat[[energy_var]])))
energy[!is.finite(energy) | energy <= 0] <- NA_real_

# Intake per 1,000 kcal for foods; food-count variable remains unchanged.
for (v in food_intake_vars) {
  food_df[[v]] <- food_df[[v]] / energy * 1000
}

zscore_numeric <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

# Standardize each outcome across the entire eligible study population.
food_df[] <- lapply(food_df, zscore_numeric)

meta <- dat[, c(group_var, model1_covars), drop = FALSE]
meta[[group_var]] <- as.character(meta[[group_var]])
for (v in continuous_covars) {
  meta[[v]] <- suppressWarnings(as.numeric(as.character(meta[[v]])))
}
for (v in categorical_covars) {
  meta[[v]] <- factor(as.character(meta[[v]]))
}

# 3. Fit one food-specific linear model -----------------------

missing_result <- function(n, status) {
  list(Beta = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
       p_value = NA_real_, n_used = n, status = status)
}

fit_single_linear <- function(y, x, covariate_df = NULL) {
  model_df <- data.frame(y = y, x = x)
  if (!is.null(covariate_df)) {
    covariate_names <- paste0("COV_", seq_len(ncol(covariate_df)))
    names(covariate_df) <- covariate_names
    model_df <- data.frame(model_df, covariate_df, check.names = FALSE)
  }

  # Complete cases are selected separately for each food and model.
  model_df <- model_df[complete.cases(model_df), , drop = FALSE]
  n <- nrow(model_df)
  if (n < 20) return(missing_result(n, "Too few complete cases"))
  if (length(unique(model_df$y)) < 2) {
    return(missing_result(n, "Food outcome has no variation"))
  }
  if (length(unique(model_df$x)) < 2) {
    return(missing_result(n, "Oral-condition contrast has no variation"))
  }

  for (v in names(model_df)) {
    if (is.factor(model_df[[v]])) model_df[[v]] <- droplevels(model_df[[v]])
  }

  # Exclude covariates with no variation after complete-case filtering.
  candidate_covars <- setdiff(names(model_df), c("y", "x"))
  valid_covars <- candidate_covars[vapply(candidate_covars, function(v) {
    z <- model_df[[v]]
    if (is.factor(z)) nlevels(z) >= 2 else length(unique(z)) >= 2
  }, logical(1))]

  formula_terms <- c("x", valid_covars)
  fit_formula <- as.formula(paste("y ~", paste(formula_terms, collapse = " + ")))
  fit <- tryCatch(suppressWarnings(lm(fit_formula, data = model_df)),
                  error = function(e) NULL)
  if (is.null(fit)) return(missing_result(n, "Model failed"))

  coefficients <- summary(fit)$coefficients
  if (!("x" %in% rownames(coefficients))) {
    return(missing_result(n, "Contrast coefficient not estimable"))
  }

  beta <- as.numeric(coefficients["x", "Estimate"])
  se <- as.numeric(coefficients["x", "Std. Error"])
  p <- as.numeric(coefficients["x", "Pr(>|t|)"])
  if (!is.finite(beta) || !is.finite(se)) {
    return(missing_result(n, "Invalid contrast estimate"))
  }

  residual_df <- df.residual(fit)
  critical <- if (is.finite(residual_df) && residual_df > 0) {
    qt(0.975, df = residual_df)
  } else {
    1.96
  }

  list(Beta = beta, CI_low = beta - critical * se,
       CI_high = beta + critical * se, p_value = p,
       n_used = n, status = "OK")
}

# 4. Run each oral-condition contrast -------------------------

run_contrast <- function(target_level, adjusted = FALSE) {
  idx <- which(meta[[group_var]] %in% c(reference_level, target_level))
  x <- as.integer(meta[[group_var]][idx] == target_level)
  covariates <- if (adjusted) meta[idx, model1_covars, drop = FALSE] else NULL

  results <- lapply(food_vars, function(food) {
    fit <- fit_single_linear(food_df[[food]][idx], x, covariates)
    data.frame(
      Food = food,
      Beta = fit$Beta,
      Beta_95CI = if (is.na(fit$Beta)) NA_character_ else sprintf(
        "%.3f (%.3f, %.3f)", fit$Beta, fit$CI_low, fit$CI_high
      ),
      p_value = fit$p_value,
      n_used = fit$n_used,
      status = fit$status,
      CI_low = fit$CI_low,
      CI_high = fit$CI_high,
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, results)
  # BH-FDR within the current contrast and model, across all food outcomes.
  result_df$FDR <- p.adjust(result_df$p_value, method = "BH")
  result_df[, c("Food", "Beta", "Beta_95CI", "p_value", "FDR",
                "n_used", "status", "CI_low", "CI_high")]
}

all_results <- list()
for (level in target_levels) {
  all_results[[paste0("Unadjusted_", level, "_vs_1")]] <-
    run_contrast(level, adjusted = FALSE)
  all_results[[paste0("Model1_", level, "_vs_1")]] <-
    run_contrast(level, adjusted = TRUE)
}

# 5. Export both models to one Excel workbook -----------------

all_results[["Model Information"]] <- data.frame(
  Setting = c("Reference", "Contrasts", "Model 1 covariates",
              "Food intake adjustment", "Food count adjustment",
              "Outcome scaling", "Multiple testing"),
  Value = c("Oral_CDC = 1 (Health)", "2 vs 1; 3 vs 1; 4 vs 1",
            paste(model1_covars, collapse = ", "),
            "Food intake per 1,000 kcal", "No energy adjustment",
            "Global z-score", "BH-FDR separately within each contrast/model"),
  stringsAsFactors = FALSE
)

openxlsx::write.xlsx(all_results, file = output_file, overwrite = TRUE)
