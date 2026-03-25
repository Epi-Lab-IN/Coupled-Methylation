
# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})

# Define input file paths 
PATH_WINDOWS <- "Escritorio/Paper_scripts/Figure2/Frog/Frog_CoM_Windows"
PATH_CHIP_SEQ <- "Escritorio/chip/H3K27me3_rana_peaks.bed"
PATH_MAPPING <- "~/Escritorio/Genomes/Assembly_frog.txt"
PATH_CHROM_SIZES <- "/home/epi-lab/Escritorio/Genomes/rana_genome.chrom.sizes"

# --- 1. Load Data and Extract Target Regions ---
frog_data <- read_delim(PATH_WINDOWS, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# Calculate upper quartile (Q3) for CA methylation
q3 <- quantile(frog_data$CA_bvalues, 0.75, na.rm = TRUE)

# Filter for high CG methylation and CA methylation above Q3 (CoM+ regions)
frog_filtered <- frog_data %>% filter(CG_bvalues > 66 & CA_bvalues > q3)

# Convert target regions to a GRanges object
frog_GR <- makeGRangesFromDataFrame(frog_filtered, keep.extra.columns = TRUE)

# --- 2. Load and Map ChIP-seq Data ---
chip_data <- read_delim(PATH_CHIP_SEQ, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
map_data <- read_delim(PATH_MAPPING, delim = "\t", show_col_types = FALSE)

# Create mapping dictionary (RefSeq <-> Scaffold)
dicc <- map_data %>%
  dplyr::select(refseq = RefSeq, scaffold = Chromosome)

# Apply mapping to ChIP-seq data and safely drop unmapped regions
chip_mapped <- chip_data %>%
  left_join(dicc, by = c("chr" = "refseq")) %>%
  mutate(chr = scaffold) %>%
  dplyr::select(-scaffold) %>%
  filter(!is.na(chr))

# Create properly mapped ChIP GRanges object
chip_GR <- GRanges(
  seqnames = chip_mapped$chr,
  ranges   = IRanges(start = chip_mapped$start, end = chip_mapped$end)
)

# Count how many of our real CoM+ regions overlap with at least one ChIP-seq peak
observed_hits <- sum(countOverlaps(frog_GR, chip_GR) > 0)

# --- 3. Generate Random Background (Memory Optimized) ---
# Note: Ensure your .chrom.sizes file actually has a header, as you specified header=TRUE
chromsizes <- read.table(PATH_CHROM_SIZES, header = TRUE)
colnames(chromsizes) <- c("chr", "size")

set.seed(42)
n_random <- length(frog_GR)

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

# Count overlaps for the random background using the SAME mapped chip_GR
expected_hits <- sum(countOverlaps(random_gr, chip_GR) > 0)

# --- 4. Statistical Analysis (Fisher's Exact Test) ---
# Construct the 2x2 contingency table
a <- observed_hits                           # Real regions WITH ChIP overlap
b <- length(frog_GR) - observed_hits         # Real regions WITHOUT ChIP overlap
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
cat("        Frog ChIP-seq Enrichment Analysis (H3K27me3)      \n")
cat("======================================================\n")
print(contingency_table)
cat("\nFisher's Exact Test Results:\n")
cat("P-value    : ", signif(fisher_res$p.value, 4), "\n")
cat("Odds Ratio : ", round(fisher_res$estimate, 2), "\n")
cat("95% CI     : [", round(fisher_res$conf.int[1], 2), "-", round(fisher_res$conf.int[2], 2), "]\n")
cat("======================================================\n")