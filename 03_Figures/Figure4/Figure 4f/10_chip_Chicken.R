# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})

# Define input file paths (Update these when publishing code)
PATH_WINDOWS <- "Escritorio/Paper_scripts/Figure2/Figura2c/Data/Chicken/Chicken_CoM_Windows"
PATH_CHIP_SEQ <- "Escritorio/chip-seq/H3K27me3_pollo_peaks.bed"
PATH_CHROM_SIZES <- "/home/epi-lab/Escritorio/Animals/pollo/GRCg7b.chrom.sizes"

# --- 1. Load Data and Extract Target Regions ---
chicken_data <- read_delim(PATH_WINDOWS, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# Calculate upper quartile (Q3) for CA methylation
q3 <- quantile(chicken_data$CA_bvalues, 0.75, na.rm = TRUE)

# Filter for high CG methylation and CA methylation above Q3
chicken_filtered <- chicken_data %>% 
  filter(CG_bvalues > 66 & CA_bvalues > q3)

# Convert target regions to a GRanges object
chicken_GR <- makeGRangesFromDataFrame(chicken_filtered, keep.extra.columns = TRUE)

# --- 2. Load ChIP-seq Data and Calculate Observed Overlaps ---
chip_data <- read_delim(PATH_CHIP_SEQ, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

chip_GR <- GRanges(
  seqnames = chip_data$chr,
  ranges   = IRanges(start = chip_data$start, end = chip_data$end)
)

# Count how many of our real CoM+ regions overlap with at least one ChIP-seq peak
observed_hits <- sum(countOverlaps(chicken_GR, chip_GR) > 0)

# --- 3. Generate Random Background (Memory Optimized) ---
# Note: Your script used header=TRUE here. Ensure your .chrom.sizes file actually has a header.
chromsizes <- read.table(PATH_CHROM_SIZES, header = TRUE)
colnames(chromsizes) <- c("chr", "size")

set.seed(42)
n_random <- length(chicken_GR)

# Sample chromosomes proportionally to their true genomic size
random_chrs <- sample(chromsizes$chr, n_random, replace = TRUE, 
                      prob = chromsizes$size / sum(chromsizes$size))

# Memory-efficient random start generation ensuring 100bp window boundaries
random_starts <- sapply(random_chrs, function(c) {
  max_start <- chromsizes$size[chromsizes$chr == c] - 100
  sample(0:((max_start - 1) %/% 100), 1) * 100 + 1
})

# Construct GRanges for the random background
random_gr <- GRanges(
  seqnames = random_chrs,
  ranges = IRanges(start = random_starts, width = 100)
)

# Count overlaps for the random background
expected_hits <- sum(countOverlaps(random_gr, chip_GR) > 0)

# --- 4. Statistical Analysis (Fisher's Exact Test) ---
# Construct the 2x2 contingency table
a <- observed_hits                           # Real regions WITH ChIP overlap
b <- length(chicken_GR) - observed_hits      # Real regions WITHOUT ChIP overlap
c <- expected_hits                           # Random regions WITH ChIP overlap
d <- length(random_gr) - expected_hits       # Random regions WITHOUT ChIP overlap

contingency_table <- matrix(
  c(a, b, 
    c, d),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    "Cohort" = c("Real CoM+", "Random"),
    "H3K27me3" = c("Overlap", "No Overlap")
  )
)

# Perform Fisher's Exact Test
fisher_res <- fisher.test(contingency_table)

# Print clean summary to the console
cat("\n======================================================\n")
cat("      Chicken ChIP-seq Enrichment Analysis (H3K27me3)     \n")
cat("======================================================\n")
print(contingency_table)
cat("\nFisher's Exact Test Results:\n")
cat("P-value    : ", signif(fisher_res$p.value, 4), "\n")
cat("Odds Ratio : ", round(fisher_res$estimate, 2), "\n")
cat("95% CI     : [", round(fisher_res$conf.int[1], 2), "-", round(fisher_res$conf.int[2], 2), "]\n")
cat("======================================================\n")