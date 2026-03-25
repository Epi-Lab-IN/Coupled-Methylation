##############################################
# Chromatin state enrichment (Odds Ratios)
# for gene body or intergenic regions
##############################################

# -----------------------------
# Load required libraries
# -----------------------------
library(biomaRt)
library(GenomicRanges)
library(genomation)
library(dplyr)
library(ggplot2)
library(readr)

# -----------------------------
# Paths and parameters
# -----------------------------
path_com          <- "~/Escritorio/Resultados/Corr_with_age/COM_pos_filtered"
path_universe_windows   <- "~/Escritorio/Resultados/Corr_with_age/universe_allgenome"
path_hmm_states   <- "~/Escritorio/Database/E073Brain_Dorsolateral_Prefrontal_Cortex primaryHMM.txt" # Obtain from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1173390897_aISDAANGAGsflKQ9uB9piMeAZYot&db=hg19&c=chr9&g=hub_1623_RoadmapConsolidatedHMM
path_gene_bed     <- "~/Escritorio/Database/hg38.bed"

# -----------------------------
# Load input data
# -----------------------------
com <- read_delim(path_com, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
target_gr   <- makeGRangesFromDataFrame(com, keep.extra.columns = TRUE)
# -----------------------------
# Annotate regions
# -----------------------------
gene_features <- readTranscriptFeatures(path_gene_bed, remove.unusual = TRUE)

# -----------------------------
# Define target and universe
# -----------------------------

universe_windows_df <- read_delim(path_universe_windows, delim = "\t", trim_ws = TRUE)
universe_gr   <- makeGRangesFromDataFrame(universe_windows_df, keep.extra.columns = TRUE)

# Annotate universe and filter by same region type
ann_universe <- annotateWithGeneParts(universe_gr, gene_features)
class_univ   <- getMembers(ann_universe)
universe_gr  <- makeGRangesFromDataFrame(
  as.data.frame(universe_gr)[(class_univ[, 2] == 1 | class_univ[, 3] == 1) & class_univ[, 1] == 0, ]
)

# -----------------------------
# Load chromatin states
# -----------------------------
chromatin_states_df <- read_delim(path_hmm_states, delim = "\t", trim_ws = TRUE)
state_order <- paste0(1:15, "_")
chromatin_state_list <- lapply(state_order, function(prefix) {
  makeGRangesFromDataFrame(
    chromatin_states_df[startsWith(chromatin_states_df$type, prefix), ],
    keep.extra.columns = TRUE
  )
})
state_names <- c(
  "Active_TSS", "Flanking_Active_TSS", "Transcription_5p_3p",
  "Strong_Transcription", "Weak_Transcription", "Genic_Enhancers",
  "Enhancers", "ZNF_Repeats", "Heterochromatin",
  "Bivalent_TSS", "Flanking_Bivalent_TSS", "Bivalent_Enhancer",
  "Repressed_PolyComb", "Weak_Repressed_PolyComb", "Quiescent/Low"
)

# -----------------------------
# Compute Odds Ratios
# -----------------------------
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
    
    contingency <- matrix(c(a, b, c, d), nrow = 2)
    fisher <- fisher.test(contingency)
    
    data.frame(
      State = state_names[i],
      a = a, b = b, c = c, d = d,
      OddsRatio = fisher$estimate,
      p_value = fisher$p.value
    )
  }) %>% bind_rows()
  
  results$p_adj <- p.adjust(results$p_value, method = "fdr")
  results$Significant <- ifelse(results$p_adj < 0.05, "Yes", "No")
  results
}

results_or <- calc_or_by_state(target_gr, universe_gr, chromatin_state_list, state_names)

# -----------------------------
# Plot
# -----------------------------
plot_or <- function(results_df, title) {
  ggplot(results_df, aes(x = OddsRatio, y = reorder(State, OddsRatio), fill = Significant)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("Yes" = "steelblue", "No" = "grey70")) +
    theme_classic() +
    labs(
      title = title,
      x = "Odds Ratio (OR)",
      y = NULL,
      fill = "Significant (FDR < 0.05)"
    ) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
}

plot_or(results_or)

# -----------------------------
# Save results
# -----------------------------
out_file <- paste0("COM_pos_close10k_", region_type, "_OR.txt")
write.table(results_or, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
