# Load required libraries
library(readr)
library(GenomicRanges)
library(dplyr)

ages <- c(57, 45, 44, 43, 40, 35)

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
ca_windows <- read_delim("/home/epi-lab/Escritorio/Animals/CHH/chimpance/chimpance_ca_comunes_Mvalues.tsv",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(ca_windows)[4:9] <- c("age35_ca", "age40_ca", "age43_ca", "age44_ca", "age45_ca", "age57_ca")
ca_windows_gr <- makeGRangesFromDataFrame(ca_windows, keep.extra.columns = TRUE)

cg_windows <- read_delim("/home/epi-lab/Escritorio/Animals/CHH/chimpance/chimpance_cg_comunes_bvalues.tsv",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(cg_windows)[4:9] <- c("age35_cg", "age40_cg", "age43_cg", "age44_cg", "age45_cg", "age57_cg")
cg_windows_gr <- makeGRangesFromDataFrame(cg_windows, keep.extra.columns = TRUE)

# Identify overlaps between CA and CG windows
overlaps <- findOverlaps(ca_windows_gr, cg_windows_gr)

# Merge overlapping regions
combined_windows <- cbind(cg_windows[subjectHits(overlaps), ], ca_windows[queryHits(overlaps), ])
combined_windows$region_id <- paste(gsub("chr", "", combined_windows$chr), combined_windows$start, combined_windows$end, sep = "-")
row.names(combined_windows) <- combined_windows$region_id
combined_windows <- combined_windows[, -c(10:12, 19)]

# Filter regions with consistently high CG methylation
high_cg <- combined_windows[apply(combined_windows[, 4:9], 1, function(x) all(x > 66)), ]

# Select and normalize CA methylation values
ca_only <- dplyr::select(high_cg, age57_ca, age45_ca, age44_ca, age43_ca, age40_ca, age35_ca)
scaled_ca_only <- as.data.frame(lapply(ca_only, scale))
row.names(scaled_ca_only) <- rownames(ca_only)

# Compute correlation for high-CG windows
lm_results_filtered <- t(apply(scaled_ca_only, 1, compute_lm_metrics))
lm_results_filtered_df <- as.data.frame(lm_results_filtered)
lm_results_filtered_df <- cbind(lm_results_filtered_df, high_cg[rownames(lm_results_filtered_df), 1:3])
sig_filtered_results <- lm_results_filtered_df[lm_results_filtered_df$p.value < 0.05, ]
sig_filtered_results <- sig_filtered_results[sig_filtered_results$slope.ages >0, ]
# Save final results
write.table(sig_filtered_results, file = "CoM_pos_chimpance", sep = "\t", row.names = FALSE, quote = FALSE)

######################## CoM- Regions ###############################

# Re-read data for overlap analysis
ca_windows <- read_delim("/home/epi-lab/Escritorio/Animals/CHH/chimpance/chimpance_ca_comunes_bvalues.tsv",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(ca_windows)[4:9] <- c("age35_ca", "age40_ca", "age43_ca", "age44_ca", "age45_ca", "age57_ca")
ca_windows_gr <- makeGRangesFromDataFrame(ca_windows, keep.extra.columns = TRUE)

cg_windows <- read_delim("/home/epi-lab/Escritorio/Animals/CHH/chimpance/chimpance_cg_comunes_Mvalues.tsv",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(cg_windows)[4:9] <- c("age35_cg", "age40_cg", "age43_cg", "age44_cg", "age45_cg", "age57_cg")
cg_windows_gr <- makeGRangesFromDataFrame(cg_windows, keep.extra.columns = TRUE)

# Identify overlaps between CA and CG windows
overlaps <- findOverlaps(ca_windows_gr, cg_windows_gr)

# Merge overlapping regions
combined_windows <- cbind(cg_windows[subjectHits(overlaps), ], ca_windows[queryHits(overlaps), ])
combined_windows$region_id <- paste(gsub("chr", "", combined_windows$chr), combined_windows$start, combined_windows$end, sep = "-")
row.names(combined_windows) <- combined_windows$region_id
combined_windows <- combined_windows[, -c(10:12, 19)]
combined_windows$mean <- rowMeans(combined_windows[10:15])
quartile <- quantile(combined_windows$mean, 0.33)
# Filter regions with consistently low CA methylation
low_ca <- combined_windows[apply(combined_windows[, 10:15], 1, function(x) all(x > quartile)), ]

# Select and normalize CA methylation values
ca_only <- dplyr::select(low_ca, age57_cg, age45_cg, age44_cg, age43_cg, age40_cg, age35_cg)
scaled_ca_only <- as.data.frame(lapply(ca_only, scale))
row.names(scaled_ca_only) <- rownames(ca_only)

# Compute correlation for low CA windows
lm_results_filtered <- t(apply(scaled_ca_only, 1, compute_lm_metrics))
lm_results_filtered_df <- as.data.frame(lm_results_filtered)
lm_results_filtered_df <- cbind(lm_results_filtered_df, low_ca[rownames(lm_results_filtered_df), 1:3])
sig_filtered_results <- lm_results_filtered_df[lm_results_filtered_df$p.value < 0.05, ]
sig_filtered_results <- sig_filtered_results[sig_filtered_results$slope.ages < 0, ]
write.table(sig_filtered_results, file = "CoM_neg_chimpance", sep = "\t", row.names = FALSE, quote = FALSE)
