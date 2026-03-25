# ================================
# Libraries
# ================================
library(readr)
library(dplyr)

# ================================
# Helper functions
# ================================

# Compute difference between old and young samples
compute_diff <- function(df, old_cols, young_cols) {
  df %>%
    mutate(
      old_mean = rowMeans(across(all_of(old_cols))),
      young_mean = rowMeans(across(all_of(young_cols))),
      diff = old_mean - young_mean
    )
}

# Read and process correlation datasets
process_correlation <- function(corr_file, methylation_df, old_cols, young_cols, threshold, positive = TRUE) {
  # Merge correlation data with methylation dataframe
  df_merged <- merge(corr_file, methylation_df, by = c("chr", "start", "end"))
  
  # Compute diff if not already present
  if(!"diff" %in% colnames(df_merged)) {
    df_merged <- compute_diff(df_merged, old_cols, young_cols)
  }
  
  # Select relevant columns
  df_selected <- df_merged %>%
    select(chr, start, end, R.edad, p.value, slope.edad, diff)
  
  # Filter by threshold
  if(positive) {
    df_filtered <- df_selected %>% filter(diff > threshold)
  } else {
    df_filtered <- df_selected %>% filter(diff < -threshold)
  }
  
  return(df_filtered)
}

# ================================
# Read data
# ================================
com_pos <- read_delim("/home/epi-lab/Escritorio/Resultados/Corr_with_age/COM_pos",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)
com_neg <- read_delim("/home/epi-lab/Escritorio/Resultados/Corr_with_age/COM_neg",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)
ca <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_ca_comun_cov5_bvalues",
                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cg <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_cg_comun_cov5_bvalues",
                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# ================================
# Compute differences for methylation datasets
# ================================
old_cols <- c("y103", "y85", "y78")
young_cols <- c("y35", "y43", "y55")

ca <- compute_diff(ca, old_cols, young_cols)
cg <- compute_diff(cg, old_cols, young_cols)

# ================================
# Process positive and negative correlations
# ================================
com_pos_filtered <- process_correlation(com_pos, ca, old_cols, young_cols, threshold = 0.05, positive = TRUE)
com_neg_filtered <- process_correlation(com_neg, cg, old_cols, young_cols, threshold = 0.05, positive = FALSE)

# ================================
# Save results
# ================================
write.table(com_pos_filtered,
            file = "COM_pos_filtrated",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(com_neg_filtered,
            file = "COM_neg_filtrated",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
