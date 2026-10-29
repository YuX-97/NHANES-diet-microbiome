# ============================================================
# Cluster-specific standardized food intake differences (k = 3)
# Export the numerical matrix originally supplied to the heatmap.
# No heatmap or other visualization is generated.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
})

# 1. Input and variable names ---------------------------------

input_file <- "dietary_intake_with_kmeans_k3.xlsx"
input_sheet <- "analysis_data"
output_file <- "Food_intake_cluster_k3_standardized_differences.xlsx"

cluster_var <- "cluster_k3"
energy_var <- "Energy_intake_kcal"

# Specify variable NAMES rather than column positions.
# Number of foods reported is not divided by total energy intake.
food_intake_vars <- character(0)
food_count_var <- "number of foods reported"
food_vars <- c(food_intake_vars, food_count_var)

# 2. Read and energy-adjust dietary variables -----------------

dat <- as.data.frame(readxl::read_excel(input_file, sheet = input_sheet))
stopifnot(length(food_intake_vars) > 0)
stopifnot(all(c(cluster_var, energy_var, food_vars) %in% names(dat)))

food_raw <- as.data.frame(lapply(
  dat[, food_vars, drop = FALSE],
  function(x) suppressWarnings(as.numeric(as.character(x)))
), check.names = FALSE)
names(food_raw) <- food_vars

energy <- suppressWarnings(as.numeric(as.character(dat[[energy_var]])))
energy[!is.finite(energy) | energy <= 0] <- NA_real_

food_ea <- food_raw
for (v in food_intake_vars) {
  food_ea[[v]] <- food_raw[[v]] / (energy / 1000)
}

# 3. Overall mean and SD for each dietary variable ------------

overall_mean <- vapply(food_ea, function(x) mean(x, na.rm = TRUE), numeric(1))
overall_sd <- vapply(food_ea, function(x) sd(x, na.rm = TRUE), numeric(1))

# Variables with zero or undefined SD cannot be standardized.
valid <- is.finite(overall_sd) & overall_sd > 0
food_ea <- food_ea[, valid, drop = FALSE]
overall_mean <- overall_mean[valid]
overall_sd <- overall_sd[valid]
stopifnot(ncol(food_ea) > 0)

# 4. Cluster means and standardized differences ---------------

cluster_codes <- sort(unique(as.character(dat[[cluster_var]][
  !is.na(dat[[cluster_var]])
])))
stopifnot(length(cluster_codes) > 0)

cluster_means <- matrix(
  NA_real_,
  nrow = ncol(food_ea),
  ncol = length(cluster_codes),
  dimnames = list(names(food_ea), paste0("Cluster ", cluster_codes))
)

cluster_n <- data.frame(
  Cluster = paste0("Cluster ", cluster_codes),
  n_participants = integer(length(cluster_codes))
)

for (i in seq_along(cluster_codes)) {
  idx <- which(as.character(dat[[cluster_var]]) == cluster_codes[i])
  cluster_n$n_participants[i] <- length(idx)
  cluster_means[, i] <- vapply(
    food_ea[idx, , drop = FALSE],
    function(x) mean(x, na.rm = TRUE),
    numeric(1)
  )
}

# Heatmap input: (cluster mean - overall mean) / overall SD.
z_matrix <- sweep(cluster_means, 1, overall_mean, FUN = "-")
z_matrix <- sweep(z_matrix, 1, overall_sd, FUN = "/")

# 5. Export the complete numerical results to Excel -----------

z_scores <- data.frame(
  Food_Variable = rownames(z_matrix),
  as.data.frame(z_matrix, check.names = FALSE),
  check.names = FALSE
)

cluster_means_df <- data.frame(
  Food_Variable = rownames(cluster_means),
  as.data.frame(cluster_means, check.names = FALSE),
  check.names = FALSE
)

overall_stats <- data.frame(
  Food_Variable = names(overall_mean),
  Overall_mean = unname(overall_mean),
  Overall_SD = unname(overall_sd),
  stringsAsFactors = FALSE
)

run_settings <- data.frame(
  Setting = c("Cluster variable", "Food intake adjustment",
              "Number of foods reported", "Standardized difference"),
  Value = c(cluster_var, "Intake per 1,000 kcal",
            "No energy adjustment",
            "(Cluster mean - overall mean) / overall SD")
)

openxlsx::write.xlsx(
  list(
    Z_scores = z_scores,
    Cluster_means = cluster_means_df,
    Overall_mean_SD = overall_stats,
    Cluster_sample_sizes = cluster_n,
    Run_settings = run_settings
  ),
  file = output_file,
  overwrite = TRUE
)
