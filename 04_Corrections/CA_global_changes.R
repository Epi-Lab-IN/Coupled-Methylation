# ================================
# Load required libraries
# ================================
library(readr)
library(Hmisc)
library(ggplot2)
library(GenomicRanges)
# ================================
# Helper Functions
# ================================

# Extract age values from column names (e.g. "y35_sample" -> 35)


# ================================
# Global ca Analysis + Hoeffding
# ================================

ca <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/CHH/CA_mvalues_common.tsv",
  delim = "\t"
)
ages <- c(35,43,55,67,78,85,103)
ca_cols <- 4:10   # ajusta si tus columnas de muestras están en otras posiciones


cat("=== Global ca Analysis ===\n")
cat("N windows:", nrow(ca), "\n")
cat("Ages extracted:", ages, "\n\n")

# Mean methylation per sample across all windows
mean_per_sample <- colMeans(ca[, ca_cols], na.rm = TRUE)

cat("Mean ca methylation per sample (%):\n")
print(round(mean_per_sample , 2))
cat("\n")

# Hoeffding test: mean ca methylation ~ age
hoef <- hoeffd(mean_per_sample, ages)

cat("--- Hoeffding Test Result ---\n")
print(hoef)

# Summary table
summary_df <- data.frame(
  context        = "mca_global",
  mean_young     = mean(mean_per_sample[1:3]) ,
  mean_old       = mean(mean_per_sample[5:7]) ,
  diff_old_young = (mean(mean_per_sample[5:7]) - mean(mean_per_sample[1:3])) ,
  hoeffding_D    = hoef$D[1, 2],
  p_value        = hoef$P[1, 2]
)

cat("\n--- Summary ---\n")
print(summary_df)


plot_df <- data.frame(
  age              = ages,
  mean_methylation = mean_per_sample ,
  sample           = colnames(ca)[ca_cols]
)


ggplot(plot_df, aes(x = age, y = mean_methylation)) +
  geom_point(size = 4, color = "#3182bd") +
  geom_smooth(method = "lm", se = TRUE, color = "#e34a33", linetype = "dashed") +
  labs(
    title = "Global CA Methylation vs Age",
    x     = "Age (years)",
    y     = "Mean CA Methylation Mvalues"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text  = element_text(color = "black", size = 14),
    axis.title = element_text(size = 16),
    axis.ticks = element_line(size = 1),
    axis.line  = element_line(size = 1)
  )


############################################## Lo mismo pero considerando solo las regiones dentro de CoM+

# ================================
# Load required libraries
# ================================
library(readr)
library(Hmisc)

# ================================
# Global ca Analysis + Hoeffding
# ================================

ca <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/CHH/CA_mvalues_common.tsv",
  delim = "\t"
)
ca_GR <- makeGRangesFromDataFrame(ca, keep.extra.columns = TRUE)
CoM_neg <- read_delim("/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/COM_neg_filtrated",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)
CoM_neg_GR <- makeGRangesFromDataFrame(CoM_neg, keep.extra.columns = TRUE)


ca_GR <- ca_GR[ca_GR%over%CoM_neg_GR,]
ca <- as.data.frame(ca_GR)
ca <- ca[-c(4,5)]
ca_cols <- 4:10   # ajusta si tus columnas de muestras están en otras posiciones


cat("=== Global ca Analysis ===\n")
cat("N windows:", nrow(ca), "\n")
cat("Ages extracted:", ages, "\n\n")

# Mean methylation per sample across all windows
mean_per_sample <- colMeans(ca[, ca_cols], na.rm = TRUE)

cat("Mean ca methylation per sample (%):\n")
print(round(mean_per_sample , 2))
cat("\n")

# Hoeffding test: mean ca methylation ~ age
hoef <- hoeffd(mean_per_sample, ages)

cat("--- Hoeffding Test Result ---\n")
print(hoef)

# Summary table
summary_df <- data.frame(
  context        = "mca_global",
  mean_young     = mean(mean_per_sample[1:3]) ,
  mean_old       = mean(mean_per_sample[5:7]) ,
  diff_old_young = (mean(mean_per_sample[5:7]) - mean(mean_per_sample[1:3])) ,
  hoeffding_D    = hoef$D[1, 2],
  p_value        = hoef$P[1, 2]
)

cat("\n--- Summary ---\n")
print(summary_df)


plot_df <- data.frame(
  age              = ages,
  mean_methylation = mean_per_sample ,
  sample           = colnames(ca)[ca_cols]
)


ggplot(plot_df, aes(x = age, y = mean_methylation)) +
  geom_point(size = 4, color = "#3182bd") +
  geom_smooth(method = "lm", se = TRUE, color = "#e34a33", linetype = "dashed") +
  labs(
    title = "Global CA Methylation vs Age",
    x     = "Age (years)",
    y     = "Mean CA Methylation Mvalues"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text  = element_text(color = "black", size = 14),
    axis.title = element_text(size = 16),
    axis.ticks = element_line(size = 1),
    axis.line  = element_line(size = 1)
  )
