# =========================================================
# ANCOM-BC Differential Abundance Analysis
#
# Method:
#   - Standard ANCOM-BC using ANCOMBC::ancombc()
#   - A single model is fitted across all samples
#   - Oral_CDC is analyzed with reference level = "1" (oral health)
#
# Input:
#   - ReadCount_3254.tsv
#   - Column 1: SEQN
#   - Column 2: Oral_CDC
#   - Columns 3 onward: microbial read counts
#
# Output:
#   - Excel file with one sheet per comparison
#   - For each Oral_CDC level (2=caries, 3=periodontitis, 4=caries+periodontitis) vs reference level "1" (oral health):
#       Genus, log2FoldChange, pvalue, padj
#
# Notes:
#   - Count matrix is converted to taxa × samples format
#   - All-zero taxa are removed before model fitting
#   - log2FoldChange = lfc / log(2)
# =========================================================

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(phyloseq)
  library(readr)
  library(openxlsx)
  library(dplyr)
})

input_path <- "ReadCount_3254.tsv"
output_path <- "ANCOMBC_OralCDC_Results.xlsx"

data <- read_tsv(input_path, show_col_types = FALSE)

rownames(data) <- make.names(data[[1]], unique = TRUE)

classification_cols <- colnames(data)[2]

metadata_raw <- data[, classification_cols, drop = FALSE]
rownames(metadata_raw) <- rownames(data)

count_data <- data[, -(1:2)]
count_data <- as.data.frame(
  lapply(count_data, function(x) as.numeric(as.character(x)))
)
rownames(count_data) <- rownames(data)

run_ancombc_analysis <- function(group_vector, count_data, ref_level = "1") {
  results_list <- list()
  
  group_factor <- factor(group_vector)
  if (!(ref_level %in% levels(group_factor))) {
    stop(
      "Reference level '", ref_level, "' not found in group levels: ",
      paste(levels(group_factor), collapse = ", ")
    )
  }
  group_factor <- relevel(group_factor, ref = ref_level)
  
  taxa_all <- colnames(count_data)
  
  counts_matrix <- t(as.matrix(count_data))
  storage.mode(counts_matrix) <- "numeric"
  
  keep_taxa_idx <- rowSums(counts_matrix) > 0
  counts_matrix <- counts_matrix[keep_taxa_idx, , drop = FALSE]
  taxa_kept <- taxa_all[keep_taxa_idx]
  
  rownames(counts_matrix) <- taxa_kept
  
  otu <- phyloseq::otu_table(counts_matrix, taxa_are_rows = TRUE)
  
  sample_df <- data.frame(
    group = group_factor,
    stringsAsFactors = FALSE
  )
  rownames(sample_df) <- rownames(count_data)
  
  sam <- phyloseq::sample_data(sample_df)
  ps <- phyloseq::phyloseq(otu, sam)
  
  out <- ANCOMBC::ancombc(
    data = ps,
    formula = "group",
    p_adj_method = "BH",
    prv_cut = 0,
    lib_cut = 0,
    group = "group",
    struc_zero = TRUE,
    neg_lb = TRUE,
    tol = 1e-5,
    max_iter = 100,
    conserve = TRUE,
    alpha = 0.05,
    global = FALSE,
    n_cl = 1,
    verbose = FALSE
  )
  
  res_list <- out$res
  
  if (is.null(res_list$lfc) || is.null(res_list$p_val) || is.null(res_list$q_val)) {
    stop(
      "ANCOMBC result is missing lfc / p_val / q_val. Available elements in res: ",
      paste(names(res_list), collapse = ", ")
    )
  }
  
  lfc_mat <- as.data.frame(res_list$lfc)
  p_mat <- as.data.frame(res_list$p_val)
  q_mat <- as.data.frame(res_list$q_val)
  
  if (nrow(lfc_mat) != length(taxa_kept)) {
    warning(
      "Row number of lfc_mat (", nrow(lfc_mat),
      ") != length of taxa_kept (", length(taxa_kept),
      "). Check alignment."
    )
  } else {
    rownames(lfc_mat) <- taxa_kept
    rownames(p_mat) <- taxa_kept
    rownames(q_mat) <- taxa_kept
  }
  
  for (level in levels(group_factor)) {
    if (level == ref_level) next
    
    base_pattern <- paste0("group", level)
    col_candidates <- grep(base_pattern, colnames(lfc_mat), value = TRUE)
    
    if (length(col_candidates) == 0) {
      warning(
        "No coefficient column matching pattern '", base_pattern,
        "' found for level ", level, ". Available columns in lfc: ",
        paste(colnames(lfc_mat), collapse = ", ")
      )
      next
    }
    
    col_name <- col_candidates[1]
    taxa_names <- rownames(lfc_mat)
    
    log2FC <- lfc_mat[, col_name] / log(2)
    pvals <- p_mat[, col_name]
    padj <- q_mat[, col_name]
    
    df <- data.frame(
      Genus = taxa_names,
      log2FoldChange = log2FC,
      pvalue = pvals,
      padj = padj,
      stringsAsFactors = FALSE
    )
    
    results_list[[paste0("Group", level, "_vs_", ref_level)]] <- df
  }
  
  results_list
}

all_results <- list()

for (group_col in classification_cols) {
  message("Running ANCOM-BC for classification: ", group_col)
  group_results <- run_ancombc_analysis(metadata_raw[[group_col]], count_data)
  all_results[[group_col]] <- group_results
}

wb <- createWorkbook()

for (grp in names(all_results)) {
  for (sheet in names(all_results[[grp]])) {
    sheet_name <- substr(paste0(grp, "_", sheet), 1, 31)
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, all_results[[grp]][[sheet]])
  }
}

addWorksheet(wb, "Metadata")
metadata_out <- tibble(
  input_file = input_path,
  method = "ANCOMBC::ancombc",
  grouping_variable = classification_cols,
  reference_level = "1",
  sample_id_column = colnames(data)[1],
  read_count_columns = "3 onward",
  output_columns = "Genus, log2FoldChange, pvalue, padj"
)
writeData(wb, "Metadata", metadata_out)

saveWorkbook(wb, output_path, overwrite = TRUE)

message("ANCOM-BC results saved to: ", output_path)