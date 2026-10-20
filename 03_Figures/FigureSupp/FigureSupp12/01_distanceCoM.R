# Load required libraries
library(ggplot2)
library(dplyr)
library(methylKit)
library(GenomicRanges)
library(readr)

# Load methylation files for CG and CA contexts
cg_files <- list("/home/epi-lab/Escritorio/controles/cx/SRR6373911_1_bismark_bt2_pe.deduplicated.CX_report.txt_CG.gz")
ca_files <- list("/home/epi-lab/Escritorio/controles/cx/SRR6373911_1_bismark_bt2_pe.deduplicated.CX_report.txt_CA.gz")

# Read CG methylation data using methylKit
cpg <- methRead(
  cg_files,
  sample.id = list("ctrl1"),
  assembly = "mm10",
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
  sample.id = list("ctrl1"),
  assembly = "mm10",
  treatment = c(0),
  pipeline = "bismarkCytosineReport",
  context = "CHH",
  header = FALSE,
  mincov = 5,
  resolution = "base"
)

# Extract first sample and calculate % methylation for CA sites
sub_noncpg <- noncpg[[1]]
sub_noncpg <- cbind(
  sub_noncpg, 
  metilacion_noncpg = (sub_noncpg$numCs / sub_noncpg$coverage) * 100
)

# Define CA sites of interest: methylated (>0%) and from chromosome 1
ca <- sub_noncpg[
  sub_noncpg$metilacion_noncpg > 0 & sub_noncpg$chr == "chr1", 
]

# Create GRanges object for CA sites (1 bp windows)
gr_ca <- GRanges(
  seqnames = "chr1", 
  ranges = IRanges(start = ca$start, width = 1)
)

# Extract CG data and compute % methylation
sub_cpg <- cpg[[1]]
sub_cpg <- cbind(
  sub_cpg, 
  metilacioncpg = (sub_cpg$numCs / sub_cpg$coverage) * 100
)

# Select unmethylated CGs (<33%) from chromosome 1
cg_um <- sub_cpg[
  sub_cpg$chr == "chr1" & sub_cpg$metilacioncpg < 33, 
]

# Create GRanges object for unmethylated CGs
gr_cg_um <- GRanges(
  seqnames = "chr1", 
  ranges = IRanges(start = cg_um$start, width = 1)
)

# Find nearest CA site to each unmethylated CG
nearest_idx <- nearest(gr_cg_um, gr_ca)
distance <- start(gr_ca[nearest_idx]) - start(gr_cg_um)

# Combine data: unmethylated CG and their nearest CA site
closest_differences_um <- cbind(
  as.data.frame(cg_um), 
  as.data.frame(ca[nearest_idx, ]), 
  distance, 
  group = "UnMethylated"
)

# Select methylated CGs (>66%) from chromosome 1
cg_m <- sub_cpg[
  sub_cpg$chr == "chr1" & sub_cpg$metilacioncpg > 66, 
]

# Create GRanges object for methylated CGs
gr_cg_m <- GRanges(
  seqnames = "chr1", 
  ranges = IRanges(start = cg_m$start, width = 1)
)

# Find nearest CA site to each methylated CG
nearest_idx <- nearest(gr_cg_m, gr_ca)
distance <- start(gr_ca[nearest_idx]) - start(gr_cg_m)

# Combine data: methylated CG and their nearest CA site
closest_differences_m <- cbind(
  as.data.frame(cg_m), 
  as.data.frame(ca[nearest_idx, ]), 
  distance, 
  group = "Methylated"
)

# Downsample methylated CGs to match number of unmethylated CGs (to avoid bias)
set.seed(123)
sub_closest_differences_m <- closest_differences_m[
  sample(1:nrow(closest_differences_m), size = nrow(closest_differences_um), replace = FALSE), 
]

# Combine balanced methylated and unmethylated CG data
closest_differences <- rbind(
  sub_closest_differences_m,
  closest_differences_um
)

# Plot: 2D density of mCA vs distance to nearest mCG, grouped by CG methylation status
ggplot(closest_differences, aes(x = distance, y = metilacion_noncpg, color = group, fill = group)) +
  stat_density_2d(
    geom = "polygon",
    alpha = 0.4,
    contour = TRUE,
    bins = 15
  ) +
  xlim(-500, 500) +
  theme_classic() +
  scale_fill_manual(
    values = c(
      "UnMethylated" = "#F4A896", 
      "Methylated" = "#A1D6C4"
    )
  ) +
  scale_color_manual(
    values = c(
      "UnMethylated" = "#F4A896", 
      "Methylated" = "#A1D6C4"
    )
  ) +
  labs(
    x = "Distance to mCG (bp)",
    y = "%mCA",
    fill = "CG Methylation lvl",
    color = "CG Methylation lvl"
  )

