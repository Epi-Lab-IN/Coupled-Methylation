##############################################
# Age-correlated CG/CA windows: spatial clustering
# + chromatin state enrichment (Odds Ratios)
##############################################
#
# Pipeline:
#   1. Load age-correlated CG/CA windows and remove
#      windows already flagged as "common" elsewhere.
#   2. Build a size-matched random background of windows.
#   3. Annotate both sets with gene parts (promoter/genbody/intergenic)
#      and keep only the region type of interest.
#   4. Find windows that cluster spatially (<10kb from another window
#      of the same group) and test whether clustering differs between
#      the age-correlated set and the random background (Fisher test).
#   5. Take the clustered age-correlated windows as the "target" set
#      and test their enrichment in chromatin states (Odds Ratios),
#      against a matched genome-wide universe -- all in the same run,
#      no need to save/reload files between the two steps.
#
##############################################

## ----------------------------------------------------------------
## 0. CONFIG -- edit only this section
## ----------------------------------------------------------------

MTHYL_CORR_FILE    <- "Escritorio/Resultados/CHH/Corr_age/CA/CA_filtrated"
COMMON_REGIONS_FILE <- "Escritorio/Resultados/CHH/Corr_age/COM_all.txt"
UNIVERSE_FILE       <- "~/Escritorio/Resultados/CHH/Corr_age/all_regions.txt"
CHROMSIZES_FILE     <- "/home/epi-lab/Escritorio/Database/hg38.chrom.sizes"          # https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
GENE_BED_FILE       <- "/home/epi-lab/Escritorio/Database/hg38.bed"                  # UCSC Table Browser: hg38, NCBI RefSeq all
HMM_STATES_FILE     <- "~/Escritorio/Database/E073Brain_Dorsolateral_Prefrontal_Cortex primaryHMM.txt" # https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=hub_1623_RoadmapConsolidatedHMM

REGION_TYPE  <- "genbody"   # one of: "genbody", "promoter", "intergenic"
MAX_DIST_BP  <- 10002          # distance threshold used to call two windows "close"
RANDOM_SEED  <- 42

OUT_CLUSTERED_WINDOWS <- paste0(REGION_TYPE, "_CA_neg_close", MAX_DIST_BP / 1000, "k.txt")
OUT_OR_RESULTS        <- paste0(OUT_CLUSTERED_WINDOWS, "_OR.txt")

## ----------------------------------------------------------------
## 1. LIBRARIES
## ----------------------------------------------------------------

library(GenomicRanges)
library(genomation)
library(dplyr)
library(ggplot2)
library(readr)

## ----------------------------------------------------------------
## 2. HELPER FUNCTIONS
## ----------------------------------------------------------------

# Keeps rows of a genomation class matrix matching the requested region type
# (columns of getMembers() are, in order: promoter, exon, intron)
region_filter <- function(class_matrix, type) {
  switch(type,
         "genbody"    = (class_matrix[, 2] == 1 | class_matrix[, 3] == 1) & class_matrix[, 1] == 0,
         "promoter"   = class_matrix[, 1] == 1,
         "intergenic" = class_matrix[, 1] == 0 & class_matrix[, 2] == 0 & class_matrix[, 3] == 0,
         stop("Unknown region type: ", type)
  )
}

# Annotates a data.frame of windows with gene parts and keeps only
# rows of the requested region type
annotate_and_filter <- function(df, gene_features, type) {
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  class_mat <- getMembers(annotateWithGeneParts(gr, gene_features))
  df[region_filter(class_mat, type), ]
}

# Adds a "chr" prefix, sorts, and computes distance to the next window
# on the same chromosome
prepare_df <- function(df, group_label) {
  df %>%
    mutate(chr = paste0("chr", chr)) %>%
    arrange(chr, start) %>%
    group_by(chr) %>%
    mutate(distance_to_next = lead(start) - end) %>%
    ungroup() %>%
    mutate(group = group_label)
}

# Returns the subset of windows that lie within `max_dist` of at least
# one other window in the same data.frame (mutual proximity clusters)
get_close_windows <- function(df, max_dist) {
  gr <- GRanges(df$chr, IRanges(df$start, df$start + 99))
  hits <- findOverlaps(gr, gr, maxgap = max_dist, type = "any")
  hits <- hits[queryHits(hits) != subjectHits(hits)]   # drop self-overlaps
  close_idx <- unique(c(queryHits(hits), subjectHits(hits)))
  df[close_idx, ]
}

# Computes Fisher-test Odds Ratios of `target` vs `universe` overlap
# with each chromatin state
calc_or_by_state <- function(target, universe, chromatin_state_list, state_names) {
  outside_target <- setdiff(universe, target, ignore.strand = TRUE)
  
  results <- lapply(seq_along(chromatin_state_list), function(i) {
    state_gr <- chromatin_state_list[[i]]
    hits_a <- findOverlaps(target, state_gr, ignore.strand = TRUE)
    hits_b <- findOverlaps(outside_target, state_gr, ignore.strand = TRUE)
    
    a <- length(unique(queryHits(hits_a)))
    c <- length(target) - a
    b <- length(unique(queryHits(hits_b)))
    d <- length(outside_target) - b
    
    fisher <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    
    data.frame(State = state_names[i], a = a, b = b, c = c, d = d,
               OddsRatio = fisher$estimate, p_value = fisher$p.value)
  }) %>% bind_rows()
  
  results$p_adj <- p.adjust(results$p_value, method = "fdr")
  results$Significant <- ifelse(results$p_adj < 0.05, "Yes", "No")
  results
}

plot_or <- function(results_df, title) {
  ggplot(results_df, aes(x = OddsRatio, y = reorder(State, OddsRatio), fill = Significant)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("Yes" = "steelblue", "No" = "grey70")) +
    theme_classic() +
    labs(title = title, x = "Odds Ratio (OR)", y = NULL, fill = "Significant (FDR < 0.05)") +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
}

## ----------------------------------------------------------------
## 3. Load age-correlated CG windows (negative slope) and remove
##    windows already flagged as "common" elsewhere
## ----------------------------------------------------------------

meth_corr <- read_delim(MTHYL_CORR_FILE, delim = "\t", trim_ws = TRUE)
meth_corr <- meth_corr[meth_corr$slope < 0, ]

common_regions <- read_delim(COMMON_REGIONS_FILE, delim = "\t", trim_ws = TRUE)
meth_corr <- anti_join(meth_corr, common_regions, by = c("chr", "start", "end"))

## ----------------------------------------------------------------
## 4. Auxiliary data: chromosome sizes, gene annotation, universe
## ----------------------------------------------------------------

chromsizes <- read.table(CHROMSIZES_FILE, header = FALSE, col.names = c("chr", "size"))
chromsizes <- chromsizes[grepl("^chr[0-9XY]+$", chromsizes$chr), ]

gene_features <- readTranscriptFeatures(GENE_BED_FILE, remove.unusual = TRUE)

universe_df <- read_delim(UNIVERSE_FILE, delim = "\t", trim_ws = TRUE)

## ----------------------------------------------------------------
## 5. Random background windows, matched in number to meth_corr
## ----------------------------------------------------------------

set.seed(RANDOM_SEED)
n_random <- nrow(meth_corr)

random_chrs  <- sample(chromsizes$chr, n_random, replace = TRUE,
                       prob = chromsizes$size / sum(chromsizes$size))
starts_by_chr <- setNames(lapply(chromsizes$size, function(sz) seq(1, sz - 100, by = 100)), chromsizes$chr)
random_starts <- mapply(function(chr, n) sample(starts_by_chr[[chr]], n, replace = TRUE),
                        chr = random_chrs, n = 1)

random_df <- as.data.frame(GRanges(seqnames = random_chrs, ranges = IRanges(start = random_starts, width = 100)))
colnames(random_df)[1] <- "chr"
random_df$chr <- gsub("^chr", "", random_df$chr)

## ----------------------------------------------------------------
## 6. Annotate with gene parts and keep only the region type of interest
## ----------------------------------------------------------------

sig_windows     <- annotate_and_filter(meth_corr, gene_features, REGION_TYPE)
non_sig_windows <- annotate_and_filter(random_df, gene_features, REGION_TYPE)
universe_windows <- annotate_and_filter(universe_df, gene_features, REGION_TYPE)

sig_df     <- prepare_df(sig_windows, "CoM")
non_sig_df <- prepare_df(non_sig_windows, "No correlation")
combined_df <- bind_rows(sig_df, non_sig_df)

## ----------------------------------------------------------------
## 7. Distance-to-next-window histogram
## ----------------------------------------------------------------

ggplot(combined_df, aes(x = distance_to_next, fill = group)) +
  geom_histogram(binwidth = 100, position = "identity", alpha = 1) +
  scale_fill_manual(values = c("CoM" = "#f05a3c", "No correlation" = "#babfba")) +
  labs(x = "Distance (bp)", y = "Frequency", fill = NULL) +
  theme_linedraw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1.5, colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.ticks = element_line(linewidth = 1),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  xlim(0, 20000) +
  ylim(0, 2000)

## ----------------------------------------------------------------
## 8. Spatial clustering (<MAX_DIST_BP) and Fisher test
## ----------------------------------------------------------------

sig_df     <- sig_df     %>% mutate(id = paste0("sig_", row_number()))
non_sig_df <- non_sig_df %>% mutate(id = paste0("non_", row_number()))
df_all <- bind_rows(sig_df, non_sig_df)

close_sig_windows     <- get_close_windows(sig_df, MAX_DIST_BP)
close_non_sig_windows <- get_close_windows(non_sig_df, MAX_DIST_BP)

df_all <- df_all %>%
  mutate(close = id %in% c(close_sig_windows$id, close_non_sig_windows$id))

contingency <- table(df_all$group, df_all$close)
print(contingency)

fisher_result <- fisher.test(contingency[2:1, ], alternative = "greater")
print(fisher_result$p.value, digits = 20)
fisher_result
## ----------------------------------------------------------------
## 9. Save clustered target windows
## ----------------------------------------------------------------

#write.table(close_sig_windows, OUT_CLUSTERED_WINDOWS, sep = "\t", row.names = FALSE)

## ----------------------------------------------------------------
## 10. Chromatin state enrichment (Odds Ratios) of the clustered
##     age-correlated windows vs. the (region-type-matched) universe
##     -- uses close_sig_windows directly, no need to re-read it from disk
## ----------------------------------------------------------------

target_gr <- makeGRangesFromDataFrame(close_sig_windows, keep.extra.columns = TRUE)
seqlevels(target_gr) <- gsub("chr", "", seqlevels(target_gr))

universe_gr <- makeGRangesFromDataFrame(universe_windows, keep.extra.columns = TRUE)
seqlevels(universe_gr) <- gsub("chr", "", seqlevels(universe_gr))

chromatin_states_df <- read_delim(HMM_STATES_FILE, delim = "\t", trim_ws = TRUE)
state_order <- paste0(1:15, "_")
chromatin_state_list <- lapply(state_order, function(prefix) {
  makeGRangesFromDataFrame(chromatin_states_df[startsWith(chromatin_states_df$type, prefix), ],
                           keep.extra.columns = TRUE)
})
state_names <- c(
  "Active_TSS", "Flanking_Active_TSS", "Transcription_5p_3p",
  "Strong_Transcription", "Weak_Transcription", "Genic_Enhancers",
  "Enhancers", "ZNF_Repeats", "Heterochromatin",
  "Bivalent_TSS", "Flanking_Bivalent_TSS", "Bivalent_Enhancer",
  "Repressed_PolyComb", "Weak_Repressed_PolyComb", "Quiescent/Low"
)

results_or <- calc_or_by_state(target_gr, universe_gr, chromatin_state_list, state_names)

plot_or(results_or, title = paste("Chromatin state enrichment -", REGION_TYPE))

write.table(results_or, file = OUT_OR_RESULTS, sep = "\t", quote = FALSE, row.names = FALSE)
