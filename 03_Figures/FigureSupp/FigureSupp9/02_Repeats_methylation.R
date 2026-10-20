# ==========================================================
# Calculate average methylation and coverage per repeat class
# ==========================================================

# ---------------------------
# Load required libraries
# ---------------------------
library(UCSCRepeatMasker)  # Repeat annotations from UCSC
library(AnnotationHub)     # Access AnnotationHub datasets
library(readr)             # Read delimited files
library(GenomicRanges)     # Genomic ranges operations
library(dplyr)             # Data manipulation

# ---------------------------
# Load repeat annotations from UCSC RepeatMasker
# ---------------------------
ah <- AnnotationHub()
repeats <- ah[["AH99003"]]  # Load human repeatMasker annotation

# Split repeats into classes
repeat_classes <- split(repeats, repeats$repClass)
class_names <- names(repeat_classes)

# Define preferred order of repeat classes
preferred_order <- c("LINE", "SINE", "LTR", "DNA", "Retroposon",
                     "Satellite", "Simple_repeat", "Low_complexity",
                     "rRNA", "scRNA", "snRNA", "srpRNA", "Unknown", "RC", "RNA")
class_names <- intersect(preferred_order, class_names)
repeat_classes <- repeat_classes[class_names]

# Remove "chr" prefix from seqlevels to match methylation window data
seqlevels(repeat_classes) <- gsub("^chr", "", seqlevels(repeat_classes))

# ---------------------------
# Load methylation windows (coverage >= 5)
# ---------------------------
ca_df <- read_delim(
  "Escritorio/Resultados/ventanas_ca_comun_cov5_bvalues",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

cg_df <- read_delim(
  "Escritorio/Resultados/ventanas_cg_comun_cov5_bvalues",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

# ---------------------------
# Convert methylation data to GRanges with mean methylation
# ---------------------------
# CA context
ca_gr <- GRanges(
  seqnames = ca_df$chr,
  ranges = IRanges(start = ca_df$start, end = ca_df$end),
  mcols = ca_df[, -(1:3)]
)
ca_df$mean_meth <- rowMeans(ca_df[, -(1:3)], na.rm = TRUE)
ca_gr$mean_meth <- ca_df$mean_meth

# CG context
cg_gr <- GRanges(
  seqnames = cg_df$chr,
  ranges = IRanges(start = cg_df$start, end = cg_df$end),
  mcols = cg_df[, -(1:3)]
)
cg_df$mean_meth <- rowMeans(cg_df[, -(1:3)], na.rm = TRUE)
cg_gr$mean_meth <- cg_df$mean_meth

# ---------------------------
# Function: compute mean methylation per repeat class
# ---------------------------
get_mean_methylation_by_class <- function(meth_gr, repeat_classes) {
  lapply(names(repeat_classes), function(cls) {
    hits <- findOverlaps(meth_gr, repeat_classes[[cls]])
    mean_val <- mean(meth_gr$mean_meth[queryHits(hits)], na.rm = TRUE)
    data.frame(class = cls, mean_meth = mean_val)
  }) %>% bind_rows()
}

# Calculate average methylation for CA and CG contexts
ca_meth_means <- get_mean_methylation_by_class(ca_gr, repeat_classes)
cg_meth_means <- get_mean_methylation_by_class(cg_gr, repeat_classes)

# ---------------------------
# PART 2: Coverage per repeat class
# ---------------------------
# Convert methylation windows to GRanges without methylation values
ca_ranges <- GRanges(seqnames = ca_df$chr, ranges = IRanges(start = ca_df$start, end = ca_df$end))
cg_ranges <- GRanges(seqnames = cg_df$chr, ranges = IRanges(start = cg_df$start, end = cg_df$end))

# Combine CA and CG windows to avoid double-counting
all_ranges <- reduce(c(ca_ranges, cg_ranges))
total_bases <- sum(width(all_ranges))  # total bases covered across all windows

# Function: compute coverage per repeat class
get_coverage_by_class <- function(all_gr, repeat_classes, total_bases) {
  lapply(names(repeat_classes), function(cls) {
    overlaps <- findOverlaps(all_gr, repeat_classes[[cls]])
    intersected <- pintersect(all_gr[queryHits(overlaps)], repeat_classes[[cls]][subjectHits(overlaps)])
    cov <- sum(width(intersected))
    data.frame(class = cls, bases = cov, perc = 100 * cov / total_bases)
  }) %>% bind_rows()
}

# Compute coverage per repeat class
coverage_df <- get_coverage_by_class(all_ranges, repeat_classes, total_bases)

# ---------------------------
# Output results
# ---------------------------
write.table(ca_meth_means, file = "ca_meth_means", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cg_meth_means, file = "cg_meth_means", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(coverage_df, file = "coverage_by_repeat_class", sep = "\t", quote = FALSE, row.names = FALSE)
