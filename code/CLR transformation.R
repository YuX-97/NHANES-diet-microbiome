# ============================================================
# Centered log-ratio (CLR) transformation of genus-level counts
# Pseudocount: 0.005
# Procedure: zero replacement in read counts -> TSS -> CLR
# ============================================================

# 1. Input and output -----------------------------------------

input_file <- "genus_read_counts.tsv"
output_file <- "genus_CLR_pseudocount_0.005.tsv"

# 2. Read genus-level count data -------------------------------

df <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# SEQN identifies participants; all other columns contain genus counts.
stopifnot("SEQN" %in% names(df))

individual_ids <- df[["SEQN"]]
genus_columns <- setdiff(names(df), "SEQN")
stopifnot(length(genus_columns) > 0)

count_data <- as.matrix(df[, genus_columns, drop = FALSE])
storage.mode(count_data) <- "numeric"

# Check input counts before transformation.
stopifnot(
  !anyNA(count_data),
  all(is.finite(count_data)),
  all(count_data >= 0),
  all(rowSums(count_data) > 0)
)

# 3. Zero replacement, TSS normalization, and CLR ------------

pseudocount <- 0.005

# Replace zero read counts only; retain all nonzero counts.
count_data[count_data == 0] <- pseudocount

# Total sum scaling (TSS) to obtain relative abundances.
relative_abundance <- sweep(
  count_data,
  1,
  rowSums(count_data),
  FUN = "/"
)

# CLR: log(relative abundance) minus the sample-wise mean log abundance.
log_abundance <- log(relative_abundance)
clr_matrix <- sweep(
  log_abundance,
  1,
  rowMeans(log_abundance),
  FUN = "-"
)

# 4. Validate and export CLR-transformed data -----------------

# CLR values should have an approximately zero mean per sample.
stopifnot(max(abs(rowMeans(clr_matrix))) < 1e-10)

clr_df <- data.frame(
  SEQN = individual_ids,
  as.data.frame(clr_matrix, check.names = FALSE),
  check.names = FALSE
)

write.table(
  clr_df,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)