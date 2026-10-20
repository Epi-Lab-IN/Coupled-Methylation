# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(data.table)
library(sva)
library(factoextra)

# Clean the environment
rm(list = ls())

# ========================================
# LOAD DATA
# ========================================
data <- fread("/home/epi-lab/Escritorio/RNA_GSE/3_datasets_updated.txt", header = TRUE)
data <- as.data.frame(data)
rownames(data) <- data$V1
data <- data[, -1]
data <- as.data.frame(t(data))

# Standardize gender labels
data$Gender <- gsub("female", "F", data$Gender)
data$Gender <- gsub("Female", "F", data$Gender)
data$Gender <- gsub("male", "M", data$Gender)
data$Gender <- gsub("Male", "M", data$Gender)

# Check sex ratios per study
table(data$Array, data$Gender)

# ========================================
# SPLIT EXPRESSION DATA AND SAMPLE INFO
# ========================================
info <- data.frame(
  Age = as.numeric(data$Age),
  Gender = data$Gender,
  Region = data$Region,
  Array = data$Array,
  row.names = rownames(data)
)

# Remove non-numeric columns (kept in info) and convert to numeric
data <- subset(data, select = -c(Age, Gender, Region, Array))
data <- mutate_all(data, function(x) as.numeric(as.character(x)))

# Sample description
table(info$Region, info$Array)
table(info$Gender, info$Array)
aggregate(info$Age ~ info$Array, FUN = quantile,
          probs = c(0, 0.25, 0.50, 0.75, 1))
hist(info$Age)

# ========================================
# LOAD CoM+ AND CoM- GENES
# ========================================
genes_pos <- read_delim("~/Escritorio/Resultados/CHH/genes/genes_COM_pos_aggregated.txt",
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                        show_col_types = FALSE)

genes_neg <- read_delim("~/Escritorio/Resultados/CHH/genes/genes_COM_neg_aggregated.txt",
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                        show_col_types = FALSE)

# ========================================
# EXAMPLE PLOT: GENE EXPRESSION VS AGE
# ========================================
x <- data
x$Age <- info$Age

ggplot(x, aes(x = Age, y = EGFR)) +
  geom_point(size = 2, color = "grey") +
  geom_smooth(color = "red", method = lm, se = TRUE, linewidth = 1) +
  stat_cor(method = "pearson") +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Median (Z-score)") +
  xlab("Age") +
  ggtitle("EGFR (GSE36192+GSE33000+GSE15222)")

# ========================================
# CoM+ / CoM- GENES (direct on samples)
# ========================================
# Check overlap between gene lists and expression data
length(intersect(colnames(data), unique(genes_pos$gene_symbol)))
length(intersect(colnames(data), unique(genes_neg$gene_symbol)))

# CoM+ genes
control_genes_pos <- subset(data, select = intersect(colnames(data), unique(genes_pos$gene_symbol)))
control_genes_pos$Age <- info$Age

# CoM- genes
control_genes_neg <- subset(data, select = intersect(colnames(data), unique(genes_neg$gene_symbol)))
control_genes_neg$Age <- info$Age

# Define age groups (young vs old)
df_analysis <- control_genes_neg %>%
  filter(Age <= 50 | Age >= 70) %>%
  mutate(
    Age_group = case_when(
      Age <= 50 ~ "<=50",
      Age >= 70 ~ ">=70"
    )
  )

# Samples per group
df_analysis %>% count(Age_group)

# Genes to analyze
gene_columns <- setdiff(colnames(df_analysis), c("Age_group", "Age"))

# ========================================
# WILCOXON TEST PER GENE
# ========================================
results <- lapply(gene_columns, function(gene) {
  
  x_young <- df_analysis[df_analysis$Age_group == "<=50", gene, drop = TRUE]
  x_old   <- df_analysis[df_analysis$Age_group == ">=70", gene, drop = TRUE]
  
  x_young <- x_young[!is.na(x_young)]
  x_old   <- x_old[!is.na(x_old)]
  
  test <- wilcox.test(x_young, x_old, exact = FALSE)
  
  data.frame(
    gene = gene,
    n_young = length(x_young),
    n_old = length(x_old),
    median_young = median(x_young),
    median_old = median(x_old),
    mean_young = mean(x_young),
    mean_old = mean(x_old),
    difference = mean(x_old) - mean(x_young),
    p_value = test$p.value
  )
})

results <- bind_rows(results)

# Multiple testing correction (Benjamini-Hochberg)
results <- results %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

# Expression direction
results <- results %>%
  mutate(
    direction = case_when(
      FDR <= 0.05 & difference > 0 ~ "Gain",
      FDR <= 0.05 & difference < 0 ~ "Loss",
      TRUE ~ "No change"
    )
  )

results %>% count(direction)

# ========================================
# GAIN VS LOSS BIAS
# ========================================
gain_loss <- results %>% filter(direction %in% c("Gain", "Loss"))

n_gain  <- sum(gain_loss$direction == "Gain")
n_loss  <- sum(gain_loss$direction == "Loss")
n_total <- n_gain + n_loss

# Binomial test: Gain vs Loss against 50:50
binom.test(x = n_gain, n = n_total, p = 0.5, alternative = "two.sided")


# Save results (the original D:/ path is a Windows path; adjust as needed)
write.table(
  results,
  file = "~/Escritorio/Resultados/CHH/genes/CoM_neg_genes_gain_loss_expression.txt",
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
