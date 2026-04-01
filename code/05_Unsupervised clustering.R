# ============================================================
# Dietary Pattern Analysis Pipeline
#
# Input:
#   - Metadata file: metadata.xlsx
#   - Sheet: Dietary_clustering
#   - Subject ID column: SEQN
#   - Dietary variables: columns 6-43
#   - Energy column: Energy intake (kcal)
#
# Workflow:
#   1) Energy adjustment of all dietary variables per 1000 kcal
#   2) Z-score standardization
#   3) PCA on standardized dietary variables
#   4) Export PCA scores and variance explained
#   5) Evaluate clustering tendency using Hopkins statistic
#   6) Evaluate k-means stability using multi-start ARI
#   7) Evaluate bootstrap stability using clusterboot Jaccard
#   8) Compare k-means and GMM for k = 2..5
#
# Notes:
#   - set.seed(42) is used for reproducibility
#   - Clustering evaluation is based on PCA axis values
#   - By default, PC1-PC16 are used for clustering when available
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(openxlsx)
  library(FNN)
  library(cluster)
  library(mclust)
  library(fpc)
})

set.seed(42)

metadata_path <- "metadata.xlsx"
metadata_sheet <- "Dietary_clustering"

output_path <- "Dietary_PCA_Clustering_Analysis.xlsx"

diet_cols_idx <- 6:43
energy_col <- "Energy intake (kcal)"
id_col <- "SEQN"

cluster_k_values <- 2:5
n_pc_for_clustering <- 16
n_kmeans_runs <- 20
clusterboot_B <- 50

raw_data <- read_excel(
  path = metadata_path,
  sheet = metadata_sheet
)

raw_data <- as.data.frame(raw_data)

if (!id_col %in% colnames(raw_data)) {
  stop("Column 'SEQN' was not found in the metadata sheet.")
}

if (!energy_col %in% colnames(raw_data)) {
  stop(paste0("Energy column '", energy_col, "' was not found in the metadata sheet."))
}

if (ncol(raw_data) < max(diet_cols_idx)) {
  stop("The metadata sheet does not contain all specified dietary columns.")
}

diet_mat_raw <- raw_data[, diet_cols_idx, drop = FALSE]
diet_mat_raw <- as.data.frame(lapply(diet_mat_raw, function(x) suppressWarnings(as.numeric(x))))

energy_vec <- suppressWarnings(as.numeric(raw_data[[energy_col]]))
energy_div <- energy_vec / 1000
energy_div[energy_div <= 0 | !is.finite(energy_div)] <- NA_real_

diet_mat_ea <- diet_mat_raw
for (j in seq_along(diet_mat_raw)) {
  diet_mat_ea[[j]] <- diet_mat_raw[[j]] / energy_div
}

diet_mat_ea <- as.data.frame(lapply(diet_mat_ea, function(x) {
  x[is.infinite(x)] <- NA_real_
  x
}))

diet_mat_z <- as.data.frame(lapply(diet_mat_ea, function(x) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- stats::sd(x, na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) {
    rep(NA_real_, length(x))
  } else {
    (x - mu) / sdv
  }
}))

colnames(diet_mat_z) <- colnames(diet_mat_ea)

complete_idx_pca <- stats::complete.cases(diet_mat_z)
diet_mat_pca <- diet_mat_z[complete_idx_pca, , drop = FALSE]

if (nrow(diet_mat_pca) < 2) {
  stop("Fewer than 2 complete samples remain for PCA after filtering.")
}

pca_res <- prcomp(diet_mat_pca, center = FALSE, scale. = FALSE)

pc_scores <- as.data.frame(pca_res$x)
pc_scores$SEQN <- raw_data[[id_col]][complete_idx_pca]
pc_scores <- pc_scores %>% select(SEQN, everything())

pc_var <- pca_res$sdev^2
pc_var_prop <- pc_var / sum(pc_var)
pc_var_cum <- cumsum(pc_var_prop)

pc_variance_df <- data.frame(
  PC = paste0("PC", seq_along(pc_var)),
  Eigenvalue = pc_var,
  Variance_Proportion = pc_var_prop,
  Cumulative_Proportion = pc_var_cum,
  stringsAsFactors = FALSE
)

available_pc_cols <- grep("^PC[0-9]+$", colnames(pc_scores), value = TRUE)
pc_cols <- available_pc_cols[seq_len(min(n_pc_for_clustering, length(available_pc_cols)))]

if (length(pc_cols) < 2) {
  stop("Fewer than 2 principal components are available for clustering.")
}

pc_mat_all <- pc_scores[, pc_cols, drop = FALSE]
pc_mat_all <- as.data.frame(lapply(pc_mat_all, function(x) suppressWarnings(as.numeric(x))))

complete_idx_cluster <- complete.cases(pc_mat_all)
pc_mat <- pc_mat_all[complete_idx_cluster, , drop = FALSE]
seqn_complete <- pc_scores$SEQN[complete_idx_cluster]

n <- nrow(pc_mat)
d <- ncol(pc_mat)

if (n < 10) {
  stop("Too few complete samples remain for clustering analysis.")
}

hopkins_stat <- function(data, m = NULL, seed = 42) {
  set.seed(seed)
  n <- nrow(data)
  d <- ncol(data)
  
  if (is.null(m)) {
    m <- floor(0.1 * n)
  }
  
  idx <- sample(seq_len(n), m)
  X <- data[idx, , drop = FALSE]
  
  mins <- apply(data, 2, min)
  maxs <- apply(data, 2, max)
  
  Y <- matrix(NA_real_, nrow = m, ncol = d)
  for (j in seq_len(d)) {
    Y[, j] <- runif(m, min = mins[j], max = maxs[j])
  }
  
  nn_X <- get.knnx(data, X, k = 2)$nn.dist[, 2]
  nn_Y <- get.knnx(data, Y, k = 1)$nn.dist[, 1]
  
  sum(nn_Y) / (sum(nn_Y) + sum(nn_X))
}

hopkins_value <- hopkins_stat(pc_mat, m = floor(0.1 * n), seed = 42)

hopkins_interpretation <- ifelse(
  hopkins_value < 0.45,
  "No clear cluster tendency",
  ifelse(hopkins_value > 0.7, "Strong cluster tendency", "Weak to moderate cluster tendency")
)

dist_mat <- dist(pc_mat, method = "euclidean")

multi_start_results <- list()

for (k in cluster_k_values) {
  labels_mat <- matrix(NA_integer_, nrow = n, ncol = n_kmeans_runs)
  
  for (i in seq_len(n_kmeans_runs)) {
    set.seed(42 + i + k * 1000)
    km <- kmeans(pc_mat, centers = k, nstart = 1, iter.max = 100)
    labels_mat[, i] <- km$cluster
  }
  
  ari_vals <- c()
  cnt <- 1
  for (i in 1:(n_kmeans_runs - 1)) {
    for (j in (i + 1):n_kmeans_runs) {
      ari_vals[cnt] <- adjustedRandIndex(labels_mat[, i], labels_mat[, j])
      cnt <- cnt + 1
    }
  }
  
  multi_start_results[[paste0("k", k)]] <- data.frame(
    k = k,
    ARI_mean = mean(ari_vals),
    ARI_sd = sd(ari_vals),
    ARI_min = min(ari_vals),
    ARI_max = max(ari_vals)
  )
}

multi_start_df <- do.call(rbind, multi_start_results)

clusterboot_results <- list()

for (k in cluster_k_values) {
  set.seed(42 + k * 100)
  cb_res <- clusterboot(
    pc_mat,
    B = clusterboot_B,
    clustermethod = kmeansCBI,
    k = k,
    seed = 42 + k * 100,
    multipleboot = FALSE,
    bscompare = TRUE
  )
  
  clusterboot_results[[paste0("k", k)]] <- data.frame(
    k = k,
    cluster = seq_along(cb_res$bootmean),
    Jaccard = cb_res$bootmean
  )
}

clusterboot_df <- do.call(rbind, clusterboot_results)

clusterboot_summary <- clusterboot_df %>%
  group_by(k) %>%
  summarise(
    Jaccard_mean = mean(Jaccard, na.rm = TRUE),
    Jaccard_min = min(Jaccard, na.rm = TRUE),
    Jaccard_max = max(Jaccard, na.rm = TRUE),
    .groups = "drop"
  )

gmm_kmeans_results <- list()
kmeans_assignments <- list()

for (k in cluster_k_values) {
  set.seed(42 + k)
  km_best <- kmeans(pc_mat, centers = k, nstart = 50, iter.max = 200)
  km_labels <- km_best$cluster
  
  sil_km <- silhouette(km_labels, dist_mat)
  sil_km_mean <- mean(sil_km[, "sil_width"])
  
  set.seed(42 + k)
  gmm <- Mclust(pc_mat, G = k, verbose = FALSE)
  gmm_labels <- gmm$classification
  gmm_bic <- gmm$bic
  
  sil_gmm <- silhouette(gmm_labels, dist_mat)
  sil_gmm_mean <- mean(sil_gmm[, "sil_width"])
  
  ari_km_gmm <- adjustedRandIndex(km_labels, gmm_labels)
  
  gmm_kmeans_results[[paste0("k", k)]] <- data.frame(
    k = k,
    Kmeans_withinss = km_best$tot.withinss,
    Kmeans_silhouette = sil_km_mean,
    GMM_BIC = gmm_bic,
    GMM_silhouette = sil_gmm_mean,
    ARI_Kmeans_vs_GMM = ari_km_gmm
  )
  
  kmeans_assignments[[paste0("k", k)]] <- data.frame(
    SEQN = seqn_complete,
    k = k,
    cluster = km_labels
  )
}

gmm_kmeans_df <- do.call(rbind, gmm_kmeans_results)
kmeans_assignments_df <- do.call(rbind, kmeans_assignments)

analysis_overview <- data.frame(
  metadata_file = metadata_path,
  metadata_sheet = metadata_sheet,
  id_column = id_col,
  energy_column = energy_col,
  dietary_columns = paste0(min(diet_cols_idx), "-", max(diet_cols_idx)),
  n_samples_for_pca = nrow(diet_mat_pca),
  n_samples_for_clustering = n,
  n_pcs_used_for_clustering = d,
  k_values_evaluated = paste(cluster_k_values, collapse = ", "),
  clusterboot_B = clusterboot_B,
  kmeans_runs = n_kmeans_runs,
  stringsAsFactors = FALSE
)

hopkins_overview <- data.frame(
  n_samples = n,
  n_dimensions = d,
  Hopkins_statistic = hopkins_value,
  Hopkins_comment = hopkins_interpretation,
  stringsAsFactors = FALSE
)

wb <- createWorkbook()

addWorksheet(wb, "Analysis_Overview")
writeData(wb, sheet = "Analysis_Overview", x = analysis_overview)

addWorksheet(wb, "Diet_EnergyAdjusted")
writeData(wb, sheet = "Diet_EnergyAdjusted", x = cbind(SEQN = raw_data[[id_col]], diet_mat_ea))

addWorksheet(wb, "Diet_Zscore")
writeData(wb, sheet = "Diet_Zscore", x = cbind(SEQN = raw_data[[id_col]], diet_mat_z))

addWorksheet(wb, "PCA_Scores")
writeData(wb, sheet = "PCA_Scores", x = pc_scores)

addWorksheet(wb, "PCA_Variance")
writeData(wb, sheet = "PCA_Variance", x = pc_variance_df)

addWorksheet(wb, "Hopkins_Overview")
writeData(wb, sheet = "Hopkins_Overview", x = hopkins_overview)

addWorksheet(wb, "Kmeans_MultiStart_ARI")
writeData(wb, sheet = "Kmeans_MultiStart_ARI", x = multi_start_df)

addWorksheet(wb, "Clusterboot_Jaccard")
writeData(wb, sheet = "Clusterboot_Jaccard", x = clusterboot_df)

addWorksheet(wb, "Clusterboot_Summary")
writeData(wb, sheet = "Clusterboot_Summary", x = clusterboot_summary)

addWorksheet(wb, "GMM_vs_Kmeans")
writeData(wb, sheet = "GMM_vs_Kmeans", x = gmm_kmeans_df)

addWorksheet(wb, "Kmeans_Assignments")
writeData(wb, sheet = "Kmeans_Assignments", x = kmeans_assignments_df)

saveWorkbook(wb, file = output_path, overwrite = TRUE)