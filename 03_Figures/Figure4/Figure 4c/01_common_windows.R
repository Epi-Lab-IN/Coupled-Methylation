#  libraries
library(ggplot2)
library(dplyr)
library(methylKit)
library(GenomicRanges)
library(readr)

outdir <- "/path/to/Fig4c/Results/01_species_windows"

# process each animal sequentially
# Find & replace every animal name in the script to process next animal
# "frog", "zebrafish", "platypus", "shark", "chicken"

# Load methylation files for CG and CA contexts ----
cg_files <- list("/path/to/zebrafish.CX_report.txt_CG.gz")
ca_files <- list("/path/to/zebrafish.CX_report.txt_CA.gz")

# Read CG methylation data using methylKit
cpg <- methRead(
  cg_files,
  sample.id = list("zebrafish_CG"),
  assembly = "",
  treatment = c(0),
  pipeline = "bismarkCytosineReport",
  context = "CG",
  header = FALSE,
  mincov = 5,
  resolution = "base"
)

# Read CA methylation data (non-CpG context, CHH used here to include CA sites)
noncpg <- methRead(
  ca_files,
  sample.id = list("zebrafish_CA"),
  assembly = "mm10",
  treatment = c(0),
  pipeline = "bismarkCytosineReport",
  context = "CHH",
  header = FALSE,
  mincov = 5,
  resolution = "base"
)

# Create 100bp windows for CG methylation ----
cg_windows <- tileMethylCounts(
  cpg,
  win.size = 100,     # Window size
  step.size = 100,    # Step size (non-overlapping windows)
  cov.bases = 0,      # No minimum number of covered bases required (1 CG min)
  mc.cores = 12       # Use 12 cores for parallel processing
)

# Save CG windows to file
write.table(
  cg_windows,
  file = paste0(outdir,"/01_species_windows/CG_windows_zebrafish"), 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Create 100bp windows for non-CG (e.g., CA) methylation ----
non_cg_windows <- tileMethylCounts(
  noncpg,
  win.size = 100,     # Window size
  step.size = 100,    # Step size
  cov.bases = 3,      # Require at least 3 covered bases per window
  mc.cores = 12
)

# Save non-CG windows to file
write.table(
  non_cg_windows, 
  file = paste0(outdir,"/01_species_windows/CA_windows_zebrafish"), 
  sep = "\t",   
  row.names = FALSE,
  quote = FALSE
)

# Find common windows across all samples in a dataset

files <- c(
  paste0(outdir,"/01_species_windows/CG_windows_zebrafish"),
  paste0(outdir,"/01_species_windows/CA_windows_zebrafish")
)

# Read each file and convert to GRanges
gr_list <- lapply(files, function(file) {
  df <- read_delim(
    file, 
    delim = "\t", 
    escape_double = FALSE, 
    trim_ws = TRUE
  )
  makeGRangesFromDataFrame(
    df, 
    keep.extra.columns = TRUE
  )
})

# Intersect all GRanges to keep only common windows across samples
common <- gr_list[[1]]
for (i in 2:length(gr_list)) {
  common <- common[common %over% gr_list[[i]]]
}

# Save the common windows to file
write.table(
  common, 
  file = paste0(outdir,"/01_species_windows/common_windows_CG1CA3_zebrafish"), 
  sep = "\t",   
  row.names = FALSE,
  quote = FALSE
)
