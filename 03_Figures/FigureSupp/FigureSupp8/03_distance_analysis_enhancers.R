library(GenomicRanges)  # For genomic range operations
library(dplyr)          # For data manipulation
library(ggplot2)         # For plotting
library(readr)          # For reading files
library(stringr)        # For string manipulation

# ================================
# Read data
# ================================
# Read the data file containing correlation results between mCA and mCG with age
com <- read_delim(
  "Escritorio/Resultados/CHH/Corr_age/COM_pos_filtrated",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

# Read chromosome sizes
chromsizes <- read.table("/home/epi-lab/Escritorio/Database/hg38.chrom.sizes", header=FALSE)
colnames(chromsizes) <- c("chr", "size")
chromsizes <- chromsizes[grepl("^chr[0-9XY]+$", chromsizes$chr), ]

# ================================
# Generate random windows proportional to chromosome sizes
# ================================
set.seed(42)
n_random <- nrow(com)

# 1. Select chromosomes proportionally to their size
random_chrs <- sample(
  chromsizes$chr,
  n_random,
  replace = TRUE,
  prob = chromsizes$size / sum(chromsizes$size)
)

# 2. Generate all possible starts ending in 001 for each chromosome
chrom_start_list <- setNames(lapply(chromsizes$size, function(sz) seq(1, sz - 100, by = 100)), chromsizes$chr)

# 3. Assign a random start to each chromosome in a vectorized manner
random_starts <- mapply(function(chr, n) sample(chrom_start_list[[chr]], n, replace = TRUE),
                        chr = random_chrs, n = 1)

# 4. Build GRanges object
random_gr <- GRanges(
  seqnames = random_chrs,
  ranges = IRanges(start = random_starts, width = 100)
)
non_sig_sampled <- as.data.frame(random_gr)
colnames(non_sig_sampled)[1] <- "chr"
non_sig_sampled$chr <- gsub("^chr", "", non_sig_sampled$chr)

# ================================
# Load and process enhancers
# ================================
enhancers <- read_delim("/home/epi-lab/Escritorio/Database/H3K27ac_HS_64.bed", delim = "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
enhancers$chr <- str_replace(enhancers$chr, "chr", "")
enhancers_GR <- makeGRangesFromDataFrame(enhancers)

# Filter significant windows overlapping with enhancers
exons_sig <- as.data.frame(as(com, "GRanges")[as(com, "GRanges") %over% enhancers_GR])
exons_sig <- exons_sig %>%
  dplyr::select(-width, -strand) %>% 
  dplyr::rename(chr = seqnames)  

# Filter non-significant windows overlapping with enhancers
exons_non_sig <- as.data.frame(as(non_sig_sampled, "GRanges")[as(non_sig_sampled, "GRanges") %over% enhancers_GR])
exons_non_sig <- exons_non_sig %>%
  dplyr::select(-width, -strand) %>% 
  dplyr::rename(chr = seqnames)  

# ================================
# Prepare dataframes and calculate distances
# ================================
# Prepare significant dataframe: select relevant columns, prefix chromosome, sort
sig_df <- exons_sig %>%
  mutate(chr = paste0("chr", chr)) %>%
  arrange(chr, start)

# Prepare non-significant dataframe in a similar way
non_sig_df <- exons_non_sig %>%
  mutate(chr = paste0("chr", chr))

# Calculate distance to the next window within each chromosome group for non-significant data
non_sig_sampled <- non_sig_df %>%
  arrange(chr, start) %>%  # Sort first
  group_by(chr) %>%
  mutate(distance_to_next = lead(start) - end) %>%
  ungroup()

# Calculate distance to the next window within each chromosome group for significant data
sig_df <- sig_df %>%
  group_by(chr) %>%
  mutate(distance_to_next = lead(start) - end) %>%
  ungroup()

# Add group labels for plotting
non_sig_sampled <- non_sig_sampled %>% mutate(group = "No correlation")
sig_df <- sig_df %>% mutate(group = "CoM")

# Combine both dataframes for comparative plotting
combined_df <- bind_rows(sig_df, non_sig_sampled)

# ================================
# Plot histogram comparing distances
# ================================
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

# ================================
# Helper function: Get close windows (<10kb) without duplication
# ================================
get_close_windows <- function(df, max_dist = 10000) {
  gr <- GRanges(df$chr, IRanges(df$start, df$start + 99))
  
  # Find all windows within max_dist
  hits <- findOverlaps(gr, gr, maxgap = max_dist, type = "any")
  
  # Remove self-overlaps
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  
  # Extract unique indices of close windows
  close_indices <- unique(c(queryHits(hits), subjectHits(hits)))
  
  # Return the filtered dataframe without duplicates
  df[close_indices, ]
}

# ================================
# Apply clustering function to both groups
# ================================
close_sig_windows <- get_close_windows(sig_df, max_dist = 10002)
close_non_sig_windows <- get_close_windows(non_sig_sampled, max_dist = 10002)

# ================================
# Fisher's test: proportion of close windows (<10kb)
# ================================
# Add a unique ID to each window
sig_df <- sig_df %>% mutate(id = paste0("sig_", row_number()))
non_sig_sampled <- non_sig_sampled %>% mutate(id = paste0("non_", row_number()))
df_all <- bind_rows(sig_df, non_sig_sampled)

# Detect close windows within each group
close_sig <- get_close_windows(sig_df, max_dist = 10000)
close_non_sig <- get_close_windows(non_sig_sampled, max_dist = 10000)

# Mark windows that belong to close clusters
df_all <- df_all %>%
  mutate(close = id %in% c(close_sig$id, close_non_sig$id))

# Create a 2x2 contingency table: group vs closeness
contingency <- table(df_all$group, df_all$close)
print(contingency)

# Run Fisher's Exact Test
fisher_result_exact <- fisher.test(contingency[2:1, ], alternative = "greater")
print(fisher_result_exact$p.value, digits=20)