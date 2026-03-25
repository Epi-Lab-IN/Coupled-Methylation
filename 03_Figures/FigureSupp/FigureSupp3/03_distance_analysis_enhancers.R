# ==========================================================
# Compare distances between correlated (CoM) and random windows
# overlapping enhancer regions
# ==========================================================

# ---------------------------
# Load required libraries
# ---------------------------
library(GenomicRanges)  # For genomic range operations
library(dplyr)          # For data manipulation
library(enrichplot)     # For enrichment plotting
library(ggplot2)        # For plotting
library(scales)         # For scaling axes
library(readr)          # For reading delimited files
library(stringr)        # For string manipulation

# ---------------------------
# Load CoM (correlated with age) windows
# ---------------------------
com_windows <- read_delim(
  "Escritorio/Resultados/Corr_with_age/COM_pos_filtered",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

# ---------------------------
# Load chromosome sizes
# ---------------------------
chrom_sizes <- read.table("/home/epi-lab/Escritorio/hg38.chrom.sizes", header = FALSE)
colnames(chrom_sizes) <- c("chr", "size")
chrom_sizes <- chrom_sizes[grepl("^chr[0-9XY]+$", chrom_sizes$chr), ]  # Keep standard chromosomes only

# ---------------------------
# Generate random genomic windows as background
# ---------------------------
set.seed(42)
n_random <- nrow(com_windows)  # Match number of random regions to CoM windows

# 1. Sample chromosomes proportional to their length
random_chrs <- sample(
  chrom_sizes$chr,
  n_random,
  replace = TRUE,
  prob = chrom_sizes$size / sum(chrom_sizes$size)
)

# 2. Generate all possible start positions in 100 bp increments
chrom_start_list <- setNames(
  lapply(chrom_sizes$size, function(sz) seq(1, sz - 100, by = 100)),
  chrom_sizes$chr
)

# 3. Assign one random start per sampled chromosome
random_starts <- mapply(
  function(chr, n) sample(chrom_start_list[[chr]], n, replace = TRUE),
  chr = random_chrs, n = 1
)

# 4. Build GRanges object for random windows
random_gr <- GRanges(
  seqnames = random_chrs,
  ranges = IRanges(start = random_starts, width = 100)
)

# Convert to data frame and clean chromosome names
random_df <- as.data.frame(random_gr)
colnames(random_df)[1] <- "chr"
random_df$chr <- gsub("^chr", "", random_df$chr)

# ---------------------------
# Load enhancer coordinates (H3K27ac peaks)
# ---------------------------
enhancers <- read_delim(
  "/home/epi-lab/Escritorio/H3K27ac_HS_64.bed",
  delim = "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE
)
enhancers$chr <- str_replace(enhancers$chr, "chr", "")
enhancers_GR <- makeGRangesFromDataFrame(enhancers)

# ---------------------------
# Identify CoM and random windows overlapping enhancers
# ---------------------------
sig_overlaps <- as.data.frame(as(com_windows, "GRanges")[as(com_windows, "GRanges") %over% enhancers_GR])
sig_overlaps <- sig_overlaps %>%
  dplyr::select(-width, -strand) %>%
  dplyr::rename(chr = seqnames)

random_overlaps <- as.data.frame(as(random_df, "GRanges")[as(random_df, "GRanges") %over% enhancers_GR])
random_overlaps <- random_overlaps %>%
  dplyr::select(-width, -strand) %>%
  dplyr::rename(chr = seqnames)

# ---------------------------
# Prepare data frames for analysis
# ---------------------------
sig_df <- sig_overlaps %>%
  mutate(chr = paste0("chr", chr)) %>%
  arrange(chr, start)

random_df <- random_overlaps %>%
  mutate(chr = paste0("chr", chr)) %>%
  arrange(chr, start)

# ---------------------------
# Compute distance to the next window on the same chromosome
# ---------------------------
random_df <- random_df %>%
  group_by(chr) %>%
  mutate(distance_to_next = lead(start) - end) %>%
  ungroup()

sig_df <- sig_df %>%
  group_by(chr) %>%
  mutate(distance_to_next = lead(start) - end) %>%
  ungroup()

# ---------------------------
# Combine both datasets for plotting
# ---------------------------
random_df <- random_df %>% mutate(group = "No correlation")
sig_df <- sig_df %>% mutate(group = "CoM")
combined_df <- bind_rows(sig_df, random_df)

# ---------------------------
# Plot histogram of distances between consecutive windows
# ---------------------------
ggplot(combined_df, aes(x = distance_to_next, fill = group)) +
  geom_histogram(binwidth = 100, position = "identity", alpha = 1) +
  scale_fill_manual(values = c("CoM" = "#0794e3", "No correlation" = "#babfba")) +
  labs(
    x = "Distance (bp)",
    y = "Frequency",
    fill = NULL
  ) +
  theme_linedraw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 1, colour = "black"),
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
