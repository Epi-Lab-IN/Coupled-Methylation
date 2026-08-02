# ==============================================================================
# Script: Regression Stability and Significance Analysis
# Description: Computes Bootstrap Confidence Intervals and Permutation-based 
#              p-values for the linear regression slope between Sample and Slope.
# ==============================================================================

# Load required libraries
library(readr)
library(boot)

# --- Load Data ---
# Note: Ensure the file path remains accessible or use a relative path if shared.
data_path <- "/home/epi-lab/Escritorio/Paper_scripts/FigureSupp/FigureSupp1/1c/WGBS_slopes.txt"

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

# ==============================================================================
# 1. Bootstrap Analysis (Regression R-squared Stability)
# ==============================================================================

# Define the function to calculate R-squared within bootstrap replicates
reg_boot_fn <- function(data, indices) {
  sample_data <- data[indices, ]
  fit <- lm(y ~ x, data = sample_data)
  return(summary(fit)$r.squared)  # Extract R-squared
}

# Run Bootstrap with 1,000 resamples
set.seed(123) # Ensures reproducibility
boot_results <- boot(data = analysis_data, statistic = reg_boot_fn, R = 1000)

# Display bootstrap metrics (Bias, Standard Error)
print(boot_results)

# Compute 95% Percentile Confidence Interval
boot_ci_results <- boot.ci(boot_results, type = "perc")
print(boot_ci_results)

# ==============================================================================
# 2. Permutation Test (Empirical p-value)
# ==============================================================================

# Fit the true observed linear model and extract R-squared
observed_fit <- lm(y ~ x, data = analysis_data)
observed_r2 <- summary(observed_fit)$r.squared


# Run Permutation test by shuffling the 'y' (Slope) variable 1,000 times
set.seed(123)
permutation_r2 <- replicate(1000, {
  y_perm <- sample(analysis_data$y)
  summary(lm(y_perm ~ analysis_data$x))$r.squared
})

# Compute empirical p-value
# (Proportion of permuted R-squared values >= observed R-squared)
empirical_p_value <- mean(permutation_r2 >= observed_r2)

# ==============================================================================
# 3. Print Final Results
# ==============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Observed R-squared:        ", round(observed_r2, 4), "\n", sep = "")
cat("Permutation p-value:       ", empirical_p_value, "\n", sep = "")
if (!is.null(boot_ci_results)) {
  cat("95% Bootstrap CI (perc):  [", 
      round(boot_ci_results$percent[4], 4), ", ", 
      round(boot_ci_results$percent[5], 4), "]\n", sep = "")
}

