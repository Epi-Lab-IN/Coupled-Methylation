# ==========================================================
# Repeat Class Enrichment Analysis for CoM Windows
# ==========================================================

# ---------------------------
# Load required libraries
# ---------------------------
library(biomaRt)         # Query Ensembl Biomart
library(GenomicRanges)   # Genomic range operations
library(genomation)      # Annotate genomic regions
library(dplyr)           # Data manipulation
library(clusterProfiler) # Enrichment analysis
library(org.Hs.eg.db)    # Human gene annotation
library(enrichplot)      # Enrichment plotting
library(ggplot2)         # Visualization
library(scales)
library(readr)
library(UCSCRepeatMasker)
library(AnnotationHub)

# ---------------------------
# Load correlation data for mCA/mCG with age
# ---------------------------
com_df <- read_delim(
  "Escritorio/Resultados/Corr_with_age/COM_pos_filtered",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

# ---------------------------
# Load chromosome sizes
# ---------------------------
chrom_sizes <- read.table("/home/epi-lab/Escritorio/hg38.chrom.sizes", header = FALSE)
colnames(chrom_sizes) <- c("chr", "size")
chrom_sizes <- chrom_sizes[grepl("^chr[0-9XY]+$", chrom_sizes$chr), ]

# ---------------------------
# Generate random windows as background
# ---------------------------
set.seed(42)
n_random <- nrow(com_df)

# Sample chromosomes proportionally to size
random_chrs <- sample(
  chrom_sizes$chr,
  n_random,
  replace = TRUE,
  prob = chrom_sizes$size / sum(chrom_sizes$size)
)

# Generate all possible start positions in 100bp steps
chrom_start_list <- setNames(
  lapply(chrom_sizes$size, function(sz) seq(1, sz - 100, by = 100)),
  chrom_sizes$chr
)

# Assign random start to each sampled chromosome
random_starts <- mapply(function(chr, n) sample(chrom_start_list[[chr]], n, replace = TRUE),
                        chr = random_chrs, n = 1)

# Build GRanges object for background
random_gr <- GRanges(seqnames = random_chrs, ranges = IRanges(start = random_starts, width = 100))
non_sig_df <- as.data.frame(random_gr)
colnames(non_sig_df)[1] <- "chr"
non_sig_df$chr <- gsub("^chr", "", non_sig_df$chr)

# ---------------------------
# Annotate significant and non-significant windows
# ---------------------------
# Load transcript features
gene_features <- readTranscriptFeatures(path_bed_annot, remove.unusual = TRUE)

# Annotate significant windows
exons_sig <- annotateWithGeneParts(as(com_df, "GRanges"), gene_features)
class_sig <- getMembers(exons_sig)
exons_sig <- com_df[class_sig[,2] == 0 & class_sig[,3] == 0 & class_sig[,1] == 0, ]

# Annotate non-significant windows
exons_non_sig <- annotateWithGeneParts(as(non_sig_df, "GRanges"), gene_features)
class_non_sig <- getMembers(exons_non_sig)
exons_non_sig <- non_sig_df[class_non_sig[,2] == 0 & class_non_sig[,3] == 0 & class_non_sig[,1] == 0, ]

# Prepare dataframes for plotting
sig_df <- exons_sig %>%
  mutate(chr = paste0("chr", chr)) %>%
  arrange(chr, start)

non_sig_df <- exons_non_sig %>%
  mutate(chr = paste0("chr", chr)) %>%
  arrange(chr, start)

# ---------------------------
# Calculate distance to next window within each chromosome
# ---------------------------
sig_df <- sig_df %>% group_by(chr) %>%
  mutate(distance_to_next = lead(start) - end) %>%
  ungroup() %>%
  mutate(group = "CoM")

non_sig_df <- non_sig_df %>% group_by(chr) %>%
  mutate(distance_to_next = lead(start) - end) %>%
  ungroup() %>%
  mutate(group = "No correlation")

# Combine dataframes
combined_df <- bind_rows(sig_df, non_sig_df)

# ---------------------------
# Plot histogram of distances
# ---------------------------
ggplot(combined_df, aes(x = distance_to_next, fill = group)) +
  geom_histogram(binwidth = 100, position = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("CoM" = "#cc79a7", "No correlation" = "#7c807c")) +
  labs(x = "Distance (bp)", y = "Frequency", fill = NULL) +
  theme_linedraw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 1, colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 1),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  xlim(0, 50000)

# ---------------------------
# Filter close and far significant windows
# ---------------------------
sig_df$chr <- gsub("^chr", "", sig_df$chr)
close_sig_windows <- sig_df %>% filter(distance_to_next < 5002)
far_sig_windows <- sig_df %>% filter(distance_to_next > 5002)

# ---------------------------
# Load RepeatMasker annotations
# ---------------------------
ah <- AnnotationHub()
repeats <- ah[["AH99003"]]
rep_classes <- split(repeats, repeats$repClass)
class_names <- names(rep_classes)
preferred_order <- c("LINE", "SINE", "LTR", "DNA", "Retroposon",
                     "Satellite", "Simple_repeat", "Low_complexity",
                     "rRNA", "scRNA", "snRNA", "srpRNA", "Unknown", "RC", "RNA")
class_names <- intersect(preferred_order, class_names)
rep_classes <- rep_classes[class_names]
seqlevels(rep_classes) <- gsub("^chr", "", seqlevels(rep_classes))

# ---------------------------
# Define universe of windows
# ---------------------------
ca_windows_df <- read_delim(path_ca_windows, delim = "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
universe_regions <- makeGRangesFromDataFrame(ca_windows_df, keep.extra.columns = TRUE)

ann_universe <- annotateWithGeneParts(as(universe_regions, "GRanges"), gene_features)
class_universe <- getMembers(ann_universe)

# Filter for intergenic windows in universe
universe <- universe_regions[class_universe[,2] == 0 & class_universe[,3] == 0 & class_universe[,1] == 0, ]

# ---------------------------
# Function to calculate Odds Ratio (OR) vs background per repeat class
# ---------------------------
calc_or_by_repeat <- function(target, universe, rep_list, class_names) {
  outside_target <- setdiff(universe, target, ignore.strand = TRUE)
  results_df <- data.frame()
  
  for (i in seq_along(rep_list)) {
    rep_gr <- rep_list[[i]]
    
    # Overlaps with target
    hits_a <- findOverlaps(target, rep_gr, ignore.strand = TRUE)
    a <- length(unique(queryHits(hits_a)))
    
    # Target non-overlapping
    c <- length(target) - a
    
    # Overlaps outside target
    hits_b <- findOverlaps(outside_target, rep_gr, ignore.strand = TRUE)
    b <- length(unique(queryHits(hits_b)))
    
    # Outside non-overlapping
    d <- length(outside_target) - b
    
    contingency_table <- matrix(c(a, b, c, d), nrow = 2)
    fisher_test <- fisher.test(contingency_table)
    
    results_df <- rbind(
      results_df,
      data.frame(
        RepeatClass = class_names[i],
        a = a, b = b, c = c, d = d,
        OddsRatio = fisher_test$estimate,
        p_value = fisher_test$p.value
      )
    )
  }
  
  # Adjust p-values
  results_df$p_adj <- p.adjust(results_df$p_value, method = "fdr")
  results_df$Significant <- ifelse(results_df$p_adj < 0.05, "Yes", "No")
  results_df$RepeatClass <- factor(results_df$RepeatClass, levels = class_names)
  return(results_df)
}

# ---------------------------
# Calculate OR for close and far significant windows
# ---------------------------
results_close <- calc_or_by_repeat(as(close_sig_windows, "GRanges"), universe, rep_classes, class_names)
results_far   <- calc_or_by_repeat(as(far_sig_windows, "GRanges"), universe, rep_classes, class_names)

# ---------------------------
# Function to plot OR
# ---------------------------
plot_or <- function(results_df, title) {
  ggplot(results_df, aes(x = OddsRatio, y = RepeatClass, fill = Significant)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("Yes" = "steelblue", "No" = "grey70")) +
    theme_classic() +
    labs(title = title, x = "Odds Ratio (OR)", y = NULL, fill = "Significant (FDR < 0.05)") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold"))
}

# ---------------------------
# Plot results
# ---------------------------
plot_or(results_close, "Repeat Class Enrichment - Close CoM vs Background")
plot_or(results_far, "Repeat Class Enrichment - Far CoM vs Background")

# ---------------------------
# Save results
# ---------------------------
write.table(results_close, file = "COM_pos_close_10k", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_far, file = "COM_pos_far_10k", sep = "\t", quote = FALSE, row.names =
              