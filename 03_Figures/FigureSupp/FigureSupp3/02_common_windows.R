library(methylKit)
library(GenomicRanges)
# Create 100bp windows for CG methylation
cg_windows <- tileMethylCounts(
  cpg,
  win.size = 100,     # Window size
  step.size = 100,    # Step size (non-overlapping windows)
  cov.bases = 0,      # No minimum number of covered bases required
  mc.cores = 12       # Use 12 cores for parallel processing
)

# Save CG windows to file
write.table(
  cg_windows,
  file = "CG_55y", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Create 100bp windows for non-CG (e.g., CA) methylation
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
  file = "CC_55y", 
  sep = "\t",   
  row.names = FALSE
)

# Find common windows across all samples in a dataset

files <- c(
  "path/to/sample1.txt",
  "path/to/sample2.txt",
  "path/to/sample3.txt",
  "path/to/sample4.txt",
  "path/to/sample5.txt",
  "path/to/sample6.txt"
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
  file = "CC_comunes_tko_dnmt3a", 
  sep = "\t",   
  row.names = FALSE
)