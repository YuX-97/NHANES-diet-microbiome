# ============================================================
# Food intake patterns: energy adjustment -> z-score -> PCA
# -> clustering assessment (k = 2:5) -> final k-means (k = 3)
# PC1-PC16 are used for all clustering analyses.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(FNN)
  library(cluster)
  library(mclust)
  library(fpc)
})

# 1. Input, variables, and analysis settings -------------------

input_file <- "food_intake_data.xlsx"
input_sheet <- "analysis_data"
output_file <- "Food_intake_PCA_kmeans_k3.xlsx"

id_var <- "SEQN"
energy_var <- "Energy intake (kcal)"

# Supply variable NAMES from the dataset, not column positions.
# The food-count variable is z-scored but is NOT energy-adjusted.
food_intake_vars <- character(0)
food_count_var <- "number of foods reported"
diet_vars <- c(food_intake_vars, food_count_var)

random_seed <- 42
n_pc <- 16
candidate_k <- 2:5
final_k <- 3

# 2. Energy adjustment and global z-score ---------------------

raw <- as.data.frame(readxl::read_excel(input_file, sheet = input_sheet))
stopifnot(length(food_intake_vars) > 0)
stopifnot(all(c(id_var, energy_var, diet_vars) %in% names(raw)))

diet_raw <- as.data.frame(lapply(
  raw[, diet_vars, drop = FALSE],
  function(x) suppressWarnings(as.numeric(as.character(x)))
), check.names = FALSE)
names(diet_raw) <- diet_vars

energy <- suppressWarnings(as.numeric(as.character(raw[[energy_var]])))
energy[!is.finite(energy) | energy <= 0] <- NA_real_

# Convert food intakes to intake per 1,000 kcal.
# Number of foods reported remains on its original count scale.
diet_ea <- diet_raw
for (v in food_intake_vars) {
  diet_ea[[v]] <- diet_raw[[v]] / (energy / 1000)
}

zscore <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mu) / s
}

diet_z <- as.data.frame(lapply(diet_ea, zscore), check.names = FALSE)
names(diet_z) <- diet_vars

# PCA uses complete cases, following the original procedure.
complete_idx <- complete.cases(diet_z)
diet_complete <- diet_z[complete_idx, , drop = FALSE]
stopifnot(nrow(diet_complete) >= 10, ncol(diet_complete) >= n_pc)

# 3. PCA on standardized dietary variables --------------------

set.seed(random_seed)
pca <- prcomp(diet_complete, center = FALSE, scale. = FALSE)

pca_scores <- data.frame(
  SEQN = raw[[id_var]][complete_idx],
  as.data.frame(pca$x),
  check.names = FALSE
)

pc_variance <- pca$sdev^2
pca_variance <- data.frame(
  PC = paste0("PC", seq_along(pc_variance)),
  Eigenvalue = pc_variance,
  Variance_Proportion = pc_variance / sum(pc_variance),
  Cumulative_Proportion = cumsum(pc_variance / sum(pc_variance))
)

# Retain PC1-PC16 for all subsequent clustering procedures.
pc_vars <- paste0("PC", seq_len(n_pc))
stopifnot(all(pc_vars %in% names(pca_scores)))
pc_mat <- as.matrix(pca_scores[, pc_vars, drop = FALSE])
storage.mode(pc_mat) <- "numeric"
n <- nrow(pc_mat)
d <- ncol(pc_mat)
pc_dist <- dist(pc_mat, method = "euclidean")

# 4. Hopkins statistic (cluster tendency) ---------------------

hopkins_stat <- function(x, m = floor(0.1 * nrow(x)), seed = 42) {
  set.seed(seed)
  n <- nrow(x)
  dims <- ncol(x)
  idx <- sample(seq_len(n), size = m)
  X <- x[idx, , drop = FALSE]

  mins <- apply(x, 2, min)
  maxs <- apply(x, 2, max)
  Y <- matrix(NA_real_, nrow = m, ncol = dims)
  for (j in seq_len(dims)) {
    Y[, j] <- runif(m, min = mins[j], max = maxs[j])
  }

  # Use the second neighbor for observed points to exclude self-matches.
  nn_X <- FNN::get.knnx(x, X, k = 2)$nn.dist[, 2]
  nn_Y <- FNN::get.knnx(x, Y, k = 1)$nn.dist[, 1]
  sum(nn_Y) / (sum(nn_Y) + sum(nn_X))
}

hopkins <- hopkins_stat(pc_mat, seed = random_seed)
hopkins_comment <- if (hopkins < 0.45) {
  "No clear cluster tendency (near-uniform)"
} else if (hopkins > 0.7) {
  "Strong cluster tendency"
} else {
  "Weak to moderate cluster tendency"
}
hopkins_overview <- data.frame(
  n_samples = n,
  n_dimensions = d,
  Hopkins_statistic = hopkins,
  Hopkins_comment = hopkins_comment
)

# 5. Multi-start k-means stability (ARI; k = 2-5) -------------

n_runs <- 20
multi_start <- lapply(candidate_k, function(k) {
  labels <- matrix(NA_integer_, nrow = n, ncol = n_runs)
  for (i in seq_len(n_runs)) {
    set.seed(random_seed + i + k * 1000)
    km <- kmeans(pc_mat, centers = k, nstart = 1, iter.max = 100)
    labels[, i] <- km$cluster
  }

  pairs <- combn(seq_len(n_runs), 2)
  ari <- apply(pairs, 2, function(pair) {
    mclust::adjustedRandIndex(labels[, pair[1]], labels[, pair[2]])
  })

  data.frame(
    k = k,
    ARI_mean = mean(ari),
    ARI_sd = sd(ari),
    ARI_min = min(ari),
    ARI_max = max(ari)
  )
})
multi_start_ari <- do.call(rbind, multi_start)

# 6. Bootstrap stability (Jaccard; k = 2-5) -------------------

bootstrap_results <- lapply(candidate_k, function(k) {
  set.seed(random_seed + k * 100)
  fit <- fpc::clusterboot(
    pc_mat,
    B = 50,
    clustermethod = fpc::kmeansCBI,
    k = k,
    seed = random_seed + k * 100,
    multipleboot = FALSE,
    bscompare = TRUE
  )
  data.frame(
    k = k,
    cluster = seq_along(fit$bootmean),
    Jaccard = fit$bootmean
  )
})
clusterboot_jaccard <- do.call(rbind, bootstrap_results)

clusterboot_summary <- do.call(rbind, lapply(
  split(clusterboot_jaccard, clusterboot_jaccard$k),
  function(z) data.frame(
    k = z$k[1],
    Jaccard_mean = mean(z$Jaccard, na.rm = TRUE),
    Jaccard_min = min(z$Jaccard, na.rm = TRUE),
    Jaccard_max = max(z$Jaccard, na.rm = TRUE)
  )
))

# 7. Compare k-means with Gaussian mixture models ------------

gmm_vs_kmeans <- do.call(rbind, lapply(candidate_k, function(k) {
  set.seed(random_seed + k)
  km <- kmeans(pc_mat, centers = k, nstart = 50, iter.max = 200)
  km_sil <- cluster::silhouette(km$cluster, pc_dist)

  set.seed(random_seed + k)
  gmm <- mclust::Mclust(pc_mat, G = k, verbose = FALSE)
  gmm_sil <- cluster::silhouette(gmm$classification, pc_dist)

  data.frame(
    k = k,
    Kmeans_withinss = km$tot.withinss,
    Kmeans_silhouette = mean(km_sil[, "sil_width"]),
    GMM_BIC = gmm$bic,
    GMM_silhouette = mean(gmm_sil[, "sil_width"]),
    ARI_Kmeans_vs_GMM = mclust::adjustedRandIndex(
      km$cluster, gmm$classification
    )
  )
}))

# 8. Final k-means model: k = 3 -------------------------------

# k=3 is the prespecified choice following the above assessments.
set.seed(random_seed)
km_final <- kmeans(
  pc_mat, centers = final_k, nstart = 100, iter.max = 200
)

# Include every input participant; incomplete dietary data receive NA.
cluster_k3 <- data.frame(
  SEQN = raw[[id_var]],
  cluster_k3 = NA_integer_
)
cluster_k3$cluster_k3[complete_idx] <- as.integer(km_final$cluster)

kmeans_k3_summary <- data.frame(
  k = final_k,
  n_samples_used = n,
  nstart = 100,
  iter_max = 200,
  tot_withinss = km_final$tot.withinss,
  betweenss = km_final$betweenss,
  totss = km_final$totss,
  between_over_total_var = km_final$betweenss / km_final$totss
)

# 9. Export all results to one Excel workbook -----------------

settings <- data.frame(
  Setting = c(
    "Food intake adjustment", "Number of foods reported",
    "Dietary standardization", "PCA center/scale",
    "PCs used for clustering", "Candidate k", "Selected k",
    "Random seed", "Multi-start runs", "Bootstrap iterations",
    "Final k-means nstart", "Final k-means iter.max"
  ),
  Value = c(
    "Intake per 1,000 kcal", "No energy adjustment; z-score only",
    "Global z-score", "FALSE / FALSE",
    paste0("PC1-PC", n_pc), paste(candidate_k, collapse = ", "),
    as.character(final_k), as.character(random_seed),
    as.character(n_runs), "50", "100", "200"
  )
)

openxlsx::write.xlsx(
  list(
    PCA_scores = pca_scores,
    PCA_variance = pca_variance,
    Hopkins_overview = hopkins_overview,
    Kmeans_multiStart_ARI = multi_start_ari,
    Clusterboot_Jaccard = clusterboot_jaccard,
    Clusterboot_summary = clusterboot_summary,
    GMM_vs_Kmeans = gmm_vs_kmeans,
    Cluster_k3 = cluster_k3,
    Kmeans_k3_summary = kmeans_k3_summary,
    Run_settings = settings
  ),
  file = output_file,
  overwrite = TRUE
)
