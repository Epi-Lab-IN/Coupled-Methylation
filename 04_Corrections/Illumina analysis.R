###################################################################
# Methylation of ch vs cg probes and its relationship with age
#
###################################################################

## ----------------------------------------------------------------
## 0. CONFIGURATION -- edit only this section to change project
## ----------------------------------------------------------------


WORK_DIR        <- "/home/epi-lab/Escritorio/Array/"
ARRAY_INFO_FILE <- "annotationEPICtable.txt"
BETA_MATRIX_FILE <- file.path(WORK_DIR, "GSE279509_normalisedBetaMatrix.csv.gz")
SERIES_MATRIX_FILE <- file.path(WORK_DIR, "GSE279509_series_matrix.txt")

MAX_DISTANCE_BP <- 1000   # maximum ch-cg distance to consider them "paired"
P_VALUE_CUTOFF  <- 1e-5   # significance threshold for the ch~cg regression
CELL_TYPE       <- "NeuNPos"

OUT_DIR <- file.path(WORK_DIR, "figures")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)


setwd(WORK_DIR)

## ----------------------------------------------------------------
## 1. LIBRARIES
## ----------------------------------------------------------------

library(dplyr)     # select, mutate, joins
library(data.table) # fread, fast rolling joins
library(ggplot2)
library(ggpubr)    # stat_cor, stat_regline_equation
library(tidyr)     # pivot_longer
library(purrr)     # map_dbl


## ----------------------------------------------------------------
## 2. HELPER FUNCTIONS
## ----------------------------------------------------------------

# Converts a beta-value matrix to M-values, avoiding +-Inf
beta_to_mvalues <- function(beta_mat, epsilon = 1e-6) {
  beta_mat <- apply(beta_mat, 2, as.numeric)
  beta_mat[beta_mat == 0] <- epsilon
  beta_mat[beta_mat == 1] <- 1 - epsilon
  log2(beta_mat / (1 - beta_mat))
}

# Given a raw data.table (column 1 = probe ID, remaining columns = samples),
# returns an M-value data.frame with a ProbeID column
build_mvalue_df <- function(dt) {
  probe_ids <- dt[[1]]
  beta_mat  <- as.matrix(dt[, -1])
  m_mat     <- beta_to_mvalues(beta_mat)
  data.frame(ProbeID = probe_ids, m_mat, check.names = FALSE)
}

###################################################################
# Probe selection based on array annotation
###################################################################

## ----------------------------------------------------------------
## 3. Load array annotation and keep cg / ch probes
## ----------------------------------------------------------------

array_info <- read.table(ARRAY_INFO_FILE, header = TRUE, sep = "\t")
array_info$chr <- as.numeric(gsub("chr", "", array_info$chr))  # non-numeric chromosomes (X, Y, M) become NA
array_info$pos <- as.numeric(array_info$pos)

CGs <- subset(array_info, Probe_Type == "cg")
CHs <- subset(array_info, Probe_Type == "ch")

CGs <- CGs[!is.na(CGs$chr), ]
CHs <- CHs[!is.na(CHs$chr), ]

## ----------------------------------------------------------------
## 4. Pair each ch probe with its nearest cg probe (same chromosome)
##
## ----------------------------------------------------------------

cg_dt <- as.data.table(CGs)[, .(chr, pos, Name_cg = Name, pos_cg = pos)]
ch_dt <- as.data.table(CHs)[, .(chr, pos, Name_ch = Name, pos_ch = pos)]

setkey(cg_dt, chr, pos)
setkey(ch_dt, chr, pos)

min_dist <- cg_dt[ch_dt, on = .(chr, pos), roll = "nearest"]
min_dist[, pos := NULL]
min_dist[, min_distance := abs(pos_ch - pos_cg)]

min_dist <- as.data.frame(min_dist)

# Density plot of ch-cg distances
p_dist <- ggplot(min_dist, aes(x = min_distance)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = 5000, linetype = "dashed", linewidth = 0.5, col = "blue") +
  xlim(0, 50000) +
  theme_minimal() +
  labs(x = "Closest CpG (bp)", y = "Density", title = "non-CpG vs CpG Distance")
print(p_dist)
ggsave(file.path(OUT_DIR, "distancia_ch_cg.png"), p_dist, width = 6, height = 4)

message("ch-cg within <100 bp: ", sum(min_dist$min_distance < 100, na.rm = TRUE))
message("ch-cg within <1000 bp: ", sum(min_dist$min_distance < 1000, na.rm = TRUE))

# Strip probe name suffixes (e.g. trailing "_TC21")
min_dist$Name_cg <- gsub("_.*", "", min_dist$Name_cg)
min_dist$Name_ch <- gsub("_.*", "", min_dist$Name_ch)

###################################################################
# Analysis of the GSE279509 dataset
###################################################################

## ----------------------------------------------------------------
## 5. Load beta-value matrix and split into cg / ch
## ----------------------------------------------------------------

data <- fread(BETA_MATRIX_FILE, header = TRUE, sep = ",")

data_cg_raw <- data[V1 %like% "cg"]
data_ch_raw <- data[V1 %like% "ch"]

m_data_cg <- build_mvalue_df(data_cg_raw)
m_data_ch <- build_mvalue_df(data_ch_raw)

## ----------------------------------------------------------------
## 6. Filter by distance and keep cg/ch pairs present in both matrices
## ----------------------------------------------------------------

min_dist_close <- subset(min_dist, min_distance < MAX_DISTANCE_BP)
min_dist_close$list <- row.names(min_dist_close)

min_dist_close_cg <- inner_join(min_dist_close, m_data_cg, by = c("Name_cg" = "ProbeID"))
min_dist_close_ch <- inner_join(min_dist_close, m_data_ch, by = c("Name_ch" = "ProbeID"))

# Keep only the pairs that have data in both cg and ch
common_ch <- intersect(min_dist_close_cg$Name_ch, min_dist_close_ch$Name_ch)
min_dist_close_cg <- subset(min_dist_close_cg, Name_ch %in% common_ch)
min_dist_close_ch <- subset(min_dist_close_ch, Name_ch %in% common_ch)

## ----------------------------------------------------------------
## 7. Reshape to long format (one row per probe-sample) using explicit
##    column names instead of fixed indices
## ----------------------------------------------------------------
sample_cols_cg <- setdiff(names(min_dist_close_cg), c("chr", "Name_cg", "pos_cg", "Name_ch", "pos_ch", "min_distance", "list"))
sample_cols_ch <- setdiff(names(min_dist_close_ch), c("chr", "Name_cg", "pos_cg", "Name_ch", "pos_ch", "min_distance", "list"))
cg_long <- min_dist_close_cg[, c("list", sample_cols_cg)] %>%
  pivot_longer(cols = -list, names_to = "sample", values_to = "cg_value")

ch_long <- min_dist_close_ch[, c("list", sample_cols_ch)] %>%
  pivot_longer(cols = -list, names_to = "sample", values_to = "ch_value")

merged_data <- inner_join(cg_long, ch_long, by = c("list", "sample"))

## ----------------------------------------------------------------
## 8. ch ~ cg regression per sample
## ----------------------------------------------------------------

results <- merged_data %>%
  group_by(sample) %>%
  summarise(model = list(lm(ch_value ~ cg_value)), .groups = "drop") %>%
  mutate(
    slope     = purrr::map_dbl(model, ~ coef(.x)[2]),
    r_squared = purrr::map_dbl(model, ~ summary(.x)$r.squared),
    p_value   = purrr::map_dbl(model, ~ summary(.x)$coefficients[2, 4])
  ) %>%
  select(sample, r_squared, slope, p_value)

message("Samples with p < 0.05: ", sum(results$p_value < 0.05, na.rm = TRUE))

p_hist <- ggplot(results, aes(x = p_value)) +
  geom_histogram(bins = 100, fill = "steelblue") +
  coord_cartesian(xlim = c(0, 0.2)) +
  theme_minimal() +
  labs(title = "p-value distribution (ch ~ cg per sample)")
print(p_hist)
ggsave(file.path(OUT_DIR, "histograma_pvalues.png"), p_hist, width = 6, height = 4)

p_dens_cg <- ggplot(cg_long, aes(x = cg_value)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  theme_minimal() +
  labs(x = "M-values", y = "Density", title = "M-values CpG probes")
print(p_dens_cg)
ggsave(file.path(OUT_DIR, "densidad_Mvalues_cg.png"), p_dens_cg, width = 6, height = 4)

p_dens_ch <- ggplot(ch_long, aes(x = ch_value)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  theme_minimal() +
  labs(x = "M-values", y = "Density", title = "M-values CpH probes")
print(p_dens_ch)
ggsave(file.path(OUT_DIR, "densidad_Mvalues_ch.png"), p_dens_ch, width = 6, height = 4)

## ----------------------------------------------------------------
## 9. Sample metadata (age, sex, cell type)
##
## ----------------------------------------------------------------

Sample_info <- fread(SERIES_MATRIX_FILE, header = TRUE, sep = "\t", skip = 27, nrows = 32)
Sample_info <- Sample_info[c(7, 12, 13, 14), ]
Sample_info <- Sample_info[, -1]

sample_meta <- as.data.frame(t(Sample_info))
names(sample_meta) <- c("Cells", "Age", "Sex", "Sample")

sample_meta$Age    <- as.numeric(gsub("age: ", "", sample_meta$Age))
sample_meta$Sex    <- gsub("Sex: ", "", sample_meta$Sex)
sample_meta$Sample <- gsub("basename: ", "", sample_meta$Sample)

table(sample_meta$Cells)
summary(sample_meta$Age)

results_table <- inner_join(results, sample_meta, by = c("sample" = "Sample"))
table(results_table$Cells, results_table$p_value < 0.05)
sample_meta[is.na(sample_meta$Age), ]
## ----------------------------------------------------------------
## 10. Filter by cell type and significance, plot slope ~ Age
## ----------------------------------------------------------------

results_cell <- subset(results_table, Cells == CELL_TYPE)
results_cell <- subset(results_cell, p_value < P_VALUE_CUTOFF)
# results_cell$FDR <- p.adjust(results_cell$p_value, method = "fdr")  # uncomment to correct for multiple testing

summary(results_cell$Age)
table(results_cell$Sex)

plot_coord_y <- min(results_cell$slope) + (max(results_cell$slope) - min(results_cell$slope)) * 0.95
plot_coord_x <- min(results_cell$Age, na.rm = TRUE)

p_age <- ggplot(results_cell, aes(x = Age, y = slope)) +
  geom_point(color = "#035887", size = 2, alpha = 0.2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#035887", fill = "#035887", alpha = 0.15, linewidth = 1) +
  ylab("CoM slope") +
  xlab("Age") + 
  scale_x_continuous(breaks = seq(20, 100, by = 20))+
  coord_cartesian(xlim = c(15, 100))  +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(linewidth = 1),  
    axis.ticks = element_line(linewidth = 1), 
    axis.text = element_text(size = 18, color = "black"),   
    axis.title = element_text(size = 19),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 19),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") 
  )

print(p_age)
ggsave(file.path(OUT_DIR, "slope_vs_age.png"), p_age, width = 6, height = 5)