# Load required libraries
library(readr)
library(GenomicRanges)
library(dplyr)

ages <- c(35, 43, 55, 67, 78, 85, 103)

# Function to perform linear regression and extract slope, p-value, and R
compute_lm_metrics <- function(row) {
  values <- as.numeric(row)
  model <- lm(values ~ ages)
  model_summary <- summary(model)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  p_val <- model_summary$coefficients[2, 4]
  r_sq <- model_summary$r.squared
  r_value <- sign(slope) * sqrt(r_sq)
  return(c(R = r_value, p.value = p_val, slope = slope,intercept = intercept))
}

######################## CoM+ Regions ###############################

# Re-read data for overlap analysis
ca_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/CHH/EMseq_ca_comunes_Mvalues.tsv",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(ca_windows)[4:10] <- c("age35_ca", "age43_ca", "age55_ca", "age67_ca", "age78_ca", "age85_ca", "age103_ca")
ca_windows_gr <- makeGRangesFromDataFrame(ca_windows, keep.extra.columns = TRUE)

cg_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/CHH/EMseq_cg_comunes_bvalues.tsv",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(cg_windows)[4:10] <- c("age35_cg", "age43_cg", "age55_cg", "age67_cg", "age78_cg", "age85_cg", "age103_cg")
cg_windows_gr <- makeGRangesFromDataFrame(cg_windows, keep.extra.columns = TRUE)

# Identify overlaps between CA and CG windows
overlaps <- findOverlaps(ca_windows_gr, cg_windows_gr)

# Merge overlapping regions
combined_windows <- cbind(cg_windows[subjectHits(overlaps), ], ca_windows[queryHits(overlaps), ])
combined_windows$region_id <- paste(gsub("chr", "", combined_windows$chr), combined_windows$start, combined_windows$end, sep = "-")
row.names(combined_windows) <- combined_windows$region_id
combined_windows <- combined_windows[, -c(11:13, 21)]

# Filter regions with consistently high CG methylation
high_cg <- combined_windows[apply(combined_windows[, 4:10], 1, function(x) all(x > 66)), ]

# Select and normalize CA methylation values
ca_only <- dplyr::select(high_cg, age35_ca, age43_ca, age55_ca, age67_ca, age78_ca, age85_ca, age103_ca)
scaled_ca_only <- as.data.frame(lapply(ca_only, scale))
row.names(scaled_ca_only) <- rownames(ca_only)

# Compute correlation for high-CG windows
lm_results_filtered <- t(apply(scaled_ca_only, 1, compute_lm_metrics))
lm_results_filtered_df <- as.data.frame(lm_results_filtered)
lm_results_filtered_df <- cbind(lm_results_filtered_df, high_cg[rownames(lm_results_filtered_df), 1:3])
sig_filtered_results <- lm_results_filtered_df[lm_results_filtered_df$p.value < 0.05 & lm_results_filtered_df$slope.ages>0 , ]

# Save final results
write.table(sig_filtered_results, file = "CoM_pos", sep = "\t", row.names = FALSE, quote = FALSE)
dim(sig_filtered_results)
######################## CoM- Regions ###############################

# Re-read data for overlap analysis
ca_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/CHH/EMseq_ca_comunes_bvalues.tsv",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(ca_windows)[4:10] <- c("age35_ca", "age43_ca", "age55_ca", "age67_ca", "age78_ca", "age85_ca", "age103_ca")
ca_windows_gr <- makeGRangesFromDataFrame(ca_windows, keep.extra.columns = TRUE)

cg_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/CHH/EMseq_cg_comunes_Mvalues.tsv",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(cg_windows)[4:10] <- c("age35_cg", "age43_cg", "age55_cg", "age67_cg", "age78_cg", "age85_cg", "age103_cg")
cg_windows_gr <- makeGRangesFromDataFrame(cg_windows, keep.extra.columns = TRUE)

# Identify overlaps between CA and CG windows
overlaps <- findOverlaps(ca_windows_gr, cg_windows_gr)

# Merge overlapping regions
combined_windows <- cbind(cg_windows[subjectHits(overlaps), ], ca_windows[queryHits(overlaps), ])
combined_windows$region_id <- paste(gsub("chr", "", combined_windows$chr), combined_windows$start, combined_windows$end, sep = "-")
row.names(combined_windows) <- combined_windows$region_id
combined_windows <- combined_windows[, -c(11:13, 21)]
combined_windows$mean <- rowMeans(combined_windows[11:17])
quartile <- quantile(combined_windows$mean, 0.33)
# Filter regions with consistently low CA methylation
low_ca <- combined_windows[apply(combined_windows[, 11:17], 1, function(x) all(x < quartile)), ]

# Select and normalize CG methylation values
cg_only <- dplyr::select(low_ca, age35_cg, age43_cg, age55_cg, age67_cg, age78_cg, age85_cg, age103_cg)
scaled_cg_only <- as.data.frame(lapply(cg_only, scale))
row.names(scaled_cg_only) <- rownames(cg_only)

# Compute correlation for CG windows
lm_results_filtered <- t(apply(scaled_cg_only, 1, compute_lm_metrics))
lm_results_filtered_df <- as.data.frame(lm_results_filtered)
lm_results_filtered_df <- cbind(lm_results_filtered_df, low_ca[rownames(lm_results_filtered_df), 1:3])
sig_filtered_results <- lm_results_filtered_df[lm_results_filtered_df$p.value < 0.05 & lm_results_filtered_df$slope.ages<0, ]

write.table(sig_filtered_results, file = "CoM_neg", sep = "\t", row.names = FALSE, quote = FALSE)
