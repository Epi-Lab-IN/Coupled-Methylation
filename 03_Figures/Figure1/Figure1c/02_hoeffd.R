# ================================
# Load required libraries
# ================================
library(readr)
library(GenomicRanges)
library(Hmisc)

# ================================
# Helper Functions
# ================================

# Function to read and convert data to GRanges
read_to_granges <- function(path) {
  df <- read_delim(path, delim = "\t")
  makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
}

# Function to merge overlapping regions between CA and CG contexts
merge_overlaps <- function(cg_data, ca_data) {
  overlaps <- findOverlaps(ca_data, cg_data)
  merged <- cbind(
    cg[subjectHits(overlaps), ],
    ca[queryHits(overlaps), ]
  )
  merged$region_id <- paste(
    gsub("chr", "", merged$chr),
    merged$start,
    merged$end,
    sep = "-"
  )
  row.names(merged) <- merged$region_id
  merged <- merged[, -c(11:13, 21)]  # remove redundant columns
  return(merged)
}

# Function to compute Hoeffding test for subsets
compute_hoeffding <- function(subset_matrix, ages) {
  mean_values <- colMeans(subset_matrix)
  hoeffd(mean_values, ages)
}

# Function to extract age values from column names
extract_ages <- function(col_names) {
  as.numeric(sub(".*?(\\d+)y$", "\\1", col_names))
}

# Function to perform Hoeffding tests for high/intermediate/low subsets
analyze_correlation <- function(data, cg_cols, ca_cols, threshold_type = "cg") {
  ages <- extract_ages(colnames(data)[cg_cols])
  
  if (threshold_type == "cg") {
    # CG context thresholds
    high_set <- data[apply(data[, cg_cols], 1, function(x) all(x > 0.66)), ]
    intermediate_set <- data[apply(data[, cg_cols], 1, function(x) all(x >= 0.33 & x <= 0.66)), ]
    low_set <- data[apply(data[, cg_cols], 1, function(x) all(x < 0.33)), ]
    
    # In CG context we correlate CA values with age
    list(
      high = compute_hoeffding(high_set[, ca_cols], ages),
      intermediate = compute_hoeffding(intermediate_set[, ca_cols], ages),
      low = compute_hoeffding(low_set[, ca_cols], ages)
    )
    
  } else {
    # CA context thresholds based on quantiles
    mean_values <- rowMeans(data[, ca_cols])
    bot <- quantile(mean_values, 0.33)
    top <- quantile(mean_values, 0.66)
    
    high_set <- data[apply(data[, ca_cols], 1, function(x) all(x > top)), ]
    intermediate_set <- data[apply(data[, ca_cols], 1, function(x) all(x >= bot& x <= top)), ]
    low_set <- data[apply(data[, ca_cols], 1, function(x) all(x < bot)), ]
    
    # In CA context we correlate CG values with age
    list(
      high = compute_hoeffding(high_set[, cg_cols], ages),
      intermediate = compute_hoeffding(intermediate_set[, cg_cols], ages),
      low = compute_hoeffding(low_set[, cg_cols], ages)
    )
  }
}

# ================================
# mCA Context Analysis
# ================================

cg_path_bval <- "/home/epi-lab/Escritorio/Resultados/CHH/EMseq_cg_comunes_bvalues.tsv"
ca_path_mval <- "/home/epi-lab/Escritorio/Resultados/CHH/EMseq_ca_comunes_Mvalues.tsv"

cg <- read_delim(cg_path_bval, delim = "\t")
ca <- read_delim(ca_path_mval, delim = "\t")

cg_GR <- makeGRangesFromDataFrame(cg, keep.extra.columns = TRUE)
ca_GR <- makeGRangesFromDataFrame(ca, keep.extra.columns = TRUE)

cg_ca_combined <- merge_overlaps(cg_GR, ca_GR)

# Run Hoeffding correlation analysis (mCA vs mCG levels)
cg_results <- analyze_correlation(
  data = cg_ca_combined,
  cg_cols = 4:10,
  ca_cols = 11:17,
  threshold_type = "cg"
)

cg_results

# ================================
# mCG Context Analysis
# ================================

cg_path_mval <- "/home/epi-lab/Escritorio/Resultados/CHH/EMseq_cg_comunes_Mvalues.tsv"
ca_path_bval <- "/home/epi-lab/Escritorio/Resultados/CHH/EMseq_ca_comunes_bvalues.tsv"

cg <- read_delim(cg_path_mval, delim = "\t")
ca <- read_delim(ca_path_bval, delim = "\t")

cg_GR <- makeGRangesFromDataFrame(cg, keep.extra.columns = TRUE)
ca_GR <- makeGRangesFromDataFrame(ca, keep.extra.columns = TRUE)

ca_cg_combined <- merge_overlaps(cg_GR, ca_GR)

# Run Hoeffding correlation analysis (mCG vs mCA levels)
ca_results <- analyze_correlation(
  data = ca_cg_combined,
  cg_cols = 4:10,
  ca_cols = 11:17,
  threshold_type = "ca"
)

ca_results


