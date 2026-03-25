# Load required libraries
library(readr)
library(GenomicRanges)
library(dplyr)

######################## CA Correlation ###############################

# Read CA methylation data
ca_data <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_ca_comun_cov5_Mvalues",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Rename methylation columns for clarity
colnames(ca_data)[4:10] <- c("age103_ca", "age85_ca", "age75_ca", "age65_ca", "age55_ca", "age45_ca", "age35_ca")

# Create GRanges object
ca_granges <- makeGRangesFromDataFrame(ca_data, keep.extra.columns = TRUE)
# Load required libraries
library(readr)
library(GenomicRanges)
library(dplyr)

ages <- c(103, 85, 78, 67, 55, 43, 35)

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
ca_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_ca_comun_cov5_Mvalues",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(ca_windows)[4:10] <- c("age35_ca", "age43_ca", "age55_ca", "age67_ca", "age78_ca", "age85_ca", "age103_ca")
ca_windows_gr <- makeGRangesFromDataFrame(ca_windows, keep.extra.columns = TRUE)

cg_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_cg_comun_cov5_bvalues",
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
high_cg <- combined_windows[apply(combined_windows[, 4:10], 1, function(x) all(x > 0.66)), ]

# Select and normalize CA methylation values
ca_only <- dplyr::select(high_cg, age103_ca, age85_ca, age78_ca, age67_ca, age55_ca, age43_ca, age35_ca)
scaled_ca_only <- as.data.frame(lapply(ca_only, scale))
row.names(scaled_ca_only) <- rownames(ca_only)

# Compute correlation for high-CG windows
lm_results_filtered <- t(apply(scaled_ca_only, 1, compute_lm_metrics))
lm_results_filtered_df <- as.data.frame(lm_results_filtered)
lm_results_filtered_df <- cbind(lm_results_filtered_df, high_cg[rownames(lm_results_filtered_df), 1:3])
sig_filtered_results <- lm_results_filtered_df[lm_results_filtered_df$p.value < 0.05, ]

# Save final results
write.table(lm_results_filtered_df, file = "CoM_intercept", sep = "\t", row.names = FALSE, quote = FALSE)

######################## CoM- Regions ###############################

# Re-read data for overlap analysis
ca_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_ca_comun_cov5_bvalues",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(ca_windows)[4:10] <- c("age35_ca", "age43_ca", "age55_ca", "age67_ca", "age78_ca", "age85_ca", "age103_ca")
ca_windows_gr <- makeGRangesFromDataFrame(ca_windows, keep.extra.columns = TRUE)

cg_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_cg_comun_cov5_Mvalues",
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
# Filter regions with consistently high CG methylation
low_ca <- combined_windows[apply(combined_windows[, 11:17], 1, function(x) all(x > quartile)), ]

# Select and normalize CA methylation values
ca_only <- dplyr::select(low_ca, age103_cg, age85_cg, age78_cg, age67_cg, age55_cg, age43_cg, age35_cg)
scaled_ca_only <- as.data.frame(lapply(ca_only, scale))
row.names(scaled_ca_only) <- rownames(ca_only)

# Compute correlation for high-CG windows
lm_results_filtered <- t(apply(scaled_ca_only, 1, compute_lm_metrics))
lm_results_filtered_df <- as.data.frame(lm_results_filtered)
lm_results_filtered_df <- cbind(lm_results_filtered_df, high_cg[rownames(lm_results_filtered_df), 1:3])
sig_filtered_results <- lm_results_filtered_df[lm_results_filtered_df$p.value < 0.05, ]
# Generate unique region identifier
ca_data$region_id <- paste(gsub("chr", "", ca_data$chr), ca_data$start, ca_data$end, sep = "-")
row.names(ca_data) <- ca_data$region_id
ca_data$region_id <- NULL

# Normalize methylation values using Z-score
scaled_ca <- as.data.frame(lapply(ca_data[4:10], scale))
row.names(scaled_ca) <- rownames(ca_data)

# Define age vector
ages <- c(103, 85, 78, 67, 55, 43, 35)

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

# Apply linear regression function to each region
lm_results_ca <- t(apply(scaled_ca, 1, compute_lm_metrics))
lm_results_ca_df <- as.data.frame(lm_results_ca)

# Add genomic coordinates to results
lm_results_ca_df <- cbind(lm_results_ca_df, ca_data[rownames(lm_results_ca_df), 1:3])

# Filter significant results (p < 0.05)
sig_results_ca <- lm_results_ca_df[lm_results_ca_df$p.value < 0.05, ]

# Save results
write.table(lm_results_ca_df, file = "ca_intercept", sep = "\t", row.names = FALSE, quote = FALSE)


######################## CG Correlation ###############################

# Read CG methylation data
cg_data <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_cg_comun_cov5_Mvalues",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Rename methylation columns
colnames(cg_data)[4:10] <- c("age103_cg", "age85_cg", "age75_cg", "age65_cg", "age55_cg", "age45_cg", "age35_cg")

# Create GRanges object
cg_granges <- makeGRangesFromDataFrame(cg_data, keep.extra.columns = TRUE)

# Generate unique region identifier
cg_data$region_id <- paste(gsub("chr", "", cg_data$chr), cg_data$start, cg_data$end, sep = "-")
row.names(cg_data) <- cg_data$region_id
cg_data$region_id <- NULL

# Normalize methylation values
scaled_cg <- as.data.frame(lapply(cg_data[4:10], scale))
row.names(scaled_cg) <- rownames(cg_data)

# Compute linear regression
lm_results_cg <- t(apply(scaled_cg, 1, compute_lm_metrics))
lm_results_cg_df <- as.data.frame(lm_results_cg)
lm_results_cg_df <- cbind(lm_results_cg_df, cg_data[rownames(lm_results_cg_df), 1:3])
sig_results_cg <- lm_results_cg_df[lm_results_cg_df$p.value < 0.05, ]

# Save CG results
write.table(sig_results_cg, file = "cg_spearman", sep = "\t", row.names = FALSE, quote = FALSE)


######################## Filtered Shared CA-CG Regions ###############################

# Re-read data for overlap analysis
ca_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_ca_comun_cov5_Mvalues",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(ca_windows)[4:10] <- c("age35_ca", "age43_ca", "age55_ca", "age67_ca", "age78_ca", "age85_ca", "age103_ca")
ca_windows_gr <- makeGRangesFromDataFrame(ca_windows, keep.extra.columns = TRUE)

cg_windows <- read_delim("/home/epi-lab/Escritorio/Resultados/ventanas_cg_comun_cov5_Mvalues",
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
high_cg <- combined_windows[apply(combined_windows[, 4:10], 1, function(x) all(x > 0.7824086)), ]

# Select and normalize CA methylation values
ca_only <- dplyr::select(high_cg, age103_ca, age85_ca, age75_ca, age65_ca, age55_ca, age45_ca, age35_ca)
scaled_ca_only <- as.data.frame(lapply(ca_only, scale))
row.names(scaled_ca_only) <- rownames(ca_only)

# Compute correlation for high-CG windows
lm_results_filtered <- t(apply(scaled_ca_only, 1, compute_lm_metrics))
lm_results_filtered_df <- as.data.frame(lm_results_filtered)
lm_results_filtered_df <- cbind(lm_results_filtered_df, high_cg[rownames(lm_results_filtered_df), 1:3])
sig_filtered_results <- lm_results_filtered_df[lm_results_filtered_df$p.value < 0.05, ]

# Save final results
write.table(lm_results_filtered_df, file = "CoM_intercept", sep = "\t", row.names = FALSE, quote = FALSE)
