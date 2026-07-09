
# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})

# Define input file paths 
PATH_WINDOWS_DIR <- "Escritorio/Resultados/CHH/"
PATH_CHIP_SEQ <- "Escritorio/chip-seq/H3K27me3_human.bed"
PATH_CHROM_SIZES <- "/home/epi-lab/Escritorio/Database/hg38.chrom.sizes"

files <- paste0(PATH_WINDOWS_DIR, c(35, 43, 55, 67, 78, 85, 103), "y_cowindows.txt")

# --- 1. Load Data and Extract Consensus Regions ---
# Function to load samples, calculate Q3 threshold, and filter
load_and_filter <- function(file) {
  df <- read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  
  # Calculate upper quartile (Q3) for CA methylation
  q3 <- quantile(df$CA_bvalues, 0.75, na.rm = TRUE)
  
  # Filter for high CG methylation and CA methylation above Q3
  df_filtered <- df %>% 
    filter(CG_bvalues >= 66 & CA_bvalues > q3) %>%
    mutate(region_id = paste(chr, start, end, sep = "_"))
  
  return(df_filtered)
}

# Apply filtering to all 7 biological replicates
filtered_list <- lapply(files, load_and_filter)

# Intersect regions across all samples to find consensus CoM+ regions
common_regions <- Reduce(intersect, lapply(filtered_list, `[[`, "region_id"))

# Extract complete rows for consensus regions using the first sample as reference
high_cg_human <- filtered_list[[1]] %>% filter(region_id %in% common_regions)

# Format chromosomes to standard UCSC format (e.g., "chr1")
high_cg_human$chr <- paste0("chr", high_cg_human$chr)

# Convert consensus regions to a GRanges object
high_cg_human_GR <- makeGRangesFromDataFrame(high_cg_human)

# --- 2. Load ChIP-seq Data and Calculate Observed Overlaps ---
chip_data <- read_delim(PATH_CHIP_SEQ, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

chip_GR <- GRanges(
  seqnames = chip_data$chr,
  ranges   = IRanges(start = chip_data$start, end = chip_data$end)
)

# Count how many of our real CoM+ regions overlap with at least one ChIP-seq peak
observed_hits <- sum(countOverlaps(high_cg_human_GR, chip_GR) > 0)

# --- 3. Generate Random Background (Memory Optimized) ---
chromsizes <- read.table(PATH_CHROM_SIZES, header = FALSE, col.names = c("chr", "size"))

# Restrict background to standard human chromosomes
chromsizes <- chromsizes %>% filter(grepl("^chr[0-9XY]+$", chr))

set.seed(42)
n_random <- length(high_cg_human_GR)

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
b <- length(high_cg_human_GR) - observed_hits  # Real regions WITHOUT ChIP overlap
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
cat("          ChIP-seq Enrichment Analysis (H3K27me3)       \n")
cat("======================================================\n")
print(contingency_table)
cat("\nFisher's Exact Test Results:\n")
cat("P-value    : ", signif(fisher_res$p.value, 4), "\n")
cat("Odds Ratio : ", round(fisher_res$estimate, 2), "\n")
cat("95% CI     : [", round(fisher_res$conf.int[1], 2), "-", round(fisher_res$conf.int[2], 2), "]\n")
cat("======================================================\n")