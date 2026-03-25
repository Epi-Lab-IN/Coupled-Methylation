# ===============================================================
# Distance between CoM− windows and nearest CoM+ / random windows
# ===============================================================

# ---------------------------
# 0. Load required libraries
# ---------------------------
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(readr)

# ---------------------------
# 1. Load CoM windows
# ---------------------------
com_pos_df <- read_delim(
  "Escritorio/Resultados/Corr_with_age/COM_pos_filtered",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
) %>%
  mutate(chr = paste0("chr", chr))

com_neg_df <- read_delim(
  "Escritorio/Resultados/Corr_with_age/COM_neg_filtered",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
) %>%
  mutate(chr = paste0("chr", chr))

# ---------------------------
# 2. Convert data frames to GRanges
# ---------------------------
gr_com_pos <- makeGRangesFromDataFrame(com_pos_df, keep.extra.columns = TRUE)
gr_com_neg <- makeGRangesFromDataFrame(com_neg_df, keep.extra.columns = TRUE)

# ---------------------------
# 3. Compute distance from each CoM− window to the nearest CoM+ window
# ---------------------------
nearest_com_pos <- nearest(gr_com_neg, gr_com_pos)
dist_neg_to_pos <- abs(start(gr_com_pos[nearest_com_pos]) - end(gr_com_neg))

# ---------------------------
# 4. Generate random background windows
# ---------------------------
set.seed(42)

# Read chromosome sizes (hg38)
chrom_sizes <- read.table("/home/epi-lab/Escritorio/hg38.chrom.sizes", header = FALSE)
colnames(chrom_sizes) <- c("chr", "size")
chrom_sizes <- chrom_sizes[grepl("^chr[0-9XY]+$", chrom_sizes$chr), ]

# Sample random chromosomes proportional to their genomic size
random_chrs <- sample(
  chrom_sizes$chr,
  length(gr_com_neg),
  replace = TRUE,
  prob = chrom_sizes$size / sum(chrom_sizes$size)
)

# Create list of valid start positions for each chromosome
chrom_start_list <- setNames(
  lapply(chrom_sizes$size, function(sz) seq(1, sz - 100, by = 100)),
  chrom_sizes$chr
)

# Randomly sample start coordinates per chromosome
random_starts <- mapply(
  function(chr, n) sample(chrom_start_list[[chr]], n, replace = TRUE),
  chr = random_chrs, n = 1
)

# Build GRanges object with random 100 bp windows
random_gr <- GRanges(
  seqnames = random_chrs,
  ranges = IRanges(start = random_starts, width = 100)
)

# ---------------------------
# 5. Compute distance from CoM− to nearest random window
# ---------------------------
nearest_random <- nearest(gr_com_neg, random_gr)
dist_neg_to_random <- abs(start(random_gr[nearest_random]) - end(gr_com_neg))

# ---------------------------
# 6. Combine both distance distributions
# ---------------------------
df_distances <- data.frame(
  distance = c(dist_neg_to_pos, dist_neg_to_random),
  group = rep(c("CoM", "Not Correlated"),
              times = c(length(dist_neg_to_pos), length(dist_neg_to_random)))
)

# ---------------------------
# 7. Plot histogram of distances
# ---------------------------
ggplot(df_distances, aes(x = distance, fill = group)) +
  geom_histogram(binwidth = 100, position = "identity", alpha = 1) +
  scale_fill_manual(values = c("CoM" = "#f28229", "Not Correlated" = "#babfba")) +
  labs(
    x = "Distance (bp)",
    y = "Frequency",
    fill = NULL
  ) +
  theme_linedraw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 1.5, colour = "black"),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.ticks.y = element_line(size = 1),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_line(size = 1),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  xlim(0, 20000) +
  ylim(0, 300)
