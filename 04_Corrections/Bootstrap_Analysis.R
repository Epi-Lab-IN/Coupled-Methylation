# ==============================================================================
# Script: Correlation Stability and Significance Analysis
# Description: Computes Bootstrap Confidence Intervals and Permutation-based 
#              p-values for Pearson correlation between Sample and Slope.
# ==============================================================================

# Load required libraries
library(readr)
library(boot)

# --- Load Data ---
# Note: Ensure the file path remains accessible or use a relative path if shared.
data_path <- "/home/epi-lab/Escritorio/Paper_scripts/Figure1/Figura1b/Data/Emseqslopes.txt"

CA_results <- read_delim(
  file = data_path,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  show_col_types = FALSE # Suppresses verbose column parsing messages
)

# --- Prepare Data ---
# Format data directly into a clean data frame for downstream analysis
analysis_data <- data.frame(
  x = CA_results$Sample,
  y = CA_results$Slope
)

# Remove potential missing values (NA) to avoid issues in correlation calculation
analysis_data <- na.omit(analysis_data)


# ==============================================================================
# 1. Bootstrap Analysis (Correlation Stability)
# ==============================================================================

# Define the function to calculate correlation within bootstrap replicates
# Using 'method = "pearson"'; change to "spearman" if non-linear robust testing is needed.
cor_boot_fn <- function(data, indices) {
  sample_data <- data[indices, ]
  return(cor(sample_data$x, sample_data$y, method = "pearson"))
}

# Run Bootstrap with 1,000 resamples
set.seed(123) # Ensures reproducibility
boot_results <- boot(data = analysis_data, statistic = cor_boot_fn, R = 1000)

# Display bootstrap metrics (Bias, Standard Error)
print(boot_results)

# Compute 95% Percentile Confidence Interval
boot_ci_results <- boot.ci(boot_results, type = "perc")
print(boot_ci_results)


# ==============================================================================
# 2. Permutation Test (Empirical p-value)
# ==============================================================================

# Calculate the true observed Pearson correlation
observed_cor <- cor(analysis_data$x, analysis_data$y, method = "pearson")

# Run Permutation test by shuffling the 'y' (Slope) variable 1,000 times
set.seed(123)
permutation_cors <- replicate(1000, cor(analysis_data$x, sample(analysis_data$y)))

# Compute two-tailed empirical p-value
# (Proportion of permuted correlations greater than or equal to the observed one)
empirical_p_value <- mean(abs(permutation_cors) >= abs(observed_cor))


# ==============================================================================
# 3. Print Final Results
# ==============================================================================
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Observed Correlation (r): ", round(observed_cor, 4), "\n", sep = "")
cat("Permutation p-value:       ", empirical_p_value, "\n", sep = "")
if (!is.null(boot_ci_results)) {
  cat("95% Bootstrap CI (perc):  [", 
      round(boot_ci_results$percent[4], 4), ", ", 
      round(boot_ci_results$percent[5], 4), "]\n", sep = "")
}