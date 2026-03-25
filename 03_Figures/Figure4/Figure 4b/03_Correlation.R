#Load libraries
library(readr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

# Load the CG methylation windows and convert to GRanges
cpg <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/100bqW_M1_CG_min3", 
  delim = "\t", 
  escape_double = FALSE,   
  trim_ws = TRUE
)
cpg <- makeGRangesFromDataFrame(cpg, keep.extra.columns = TRUE)

# Load the common CG windows and convert to GRanges
comun_cpg <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/100bqW_M1_CG_min3", 
  delim = "\t", 
  escape_double = FALSE,   
  trim_ws = TRUE
)
cgsample <- makeGRangesFromDataFrame(comun_cpg, keep.extra.columns = TRUE)

# Keep only CG windows overlapping with the common regions
cgsample <- cpg[cpg %over% cgsample, ]
cgsample <- as.data.frame(cgsample)
cgsample$width <- NULL  # Remove unnecessary column
colnames(cgsample)[1] <- "chr"  # Rename chromosome column for consistency

# Load non-CG methylation windows
non_cg_windows <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/100bqW_M1_CA_min3",
  delim = "\t", 
  escape_double = FALSE,  
  trim_ws = TRUE
)

# Calculate % methylation and M-values for non-CG
non_cg_windows <- cbind(
  non_cg_windows, 
  CA_methylation = (non_cg_windows$numCs / non_cg_windows$coverage) * 100
)
non_cg_windows$mvalues_ca <- log2((non_cg_windows$numCs + 0.1) / (non_cg_windows$numTs + 0.1))

# Convert to GRanges
non_cg_windows_GR <- makeGRangesFromDataFrame(non_cg_windows, keep.extra.columns = TRUE)

# Load common non-CG regions and filter
comun <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/CA_comunes_10cov_min3",
  delim = "\t", 
  escape_double = FALSE,   
  trim_ws = TRUE
)
comun_GR <- makeGRangesFromDataFrame(comun, keep.extra.columns = TRUE)
non_cg_windows_GR <- non_cg_windows_GR[non_cg_windows_GR %over% comun_GR, ]
non_cg_windows <- as.data.frame(non_cg_windows_GR)
non_cg_windows$width <- NULL  # Remove unnecessary column

# Calculate CG % methylation and M-values
cgsample <- cbind(
  cgsample, 
  CG_methylation = (cgsample$numCs / cgsample$coverage) * 100
)
cgsample$mvalues_cg <- log2((cgsample$numCs + 0.1) / (cgsample$numTs + 0.1))

# Convert CG sample to GRanges
cgsample_GR <- makeGRangesFromDataFrame(cgsample, keep.extra.columns = TRUE)

# Plot density distribution of CG methylation
ggplot(cgsample, aes(x = CG_methylation)) +
  geom_density(fill = "#3182bd", color = "#1c6ca6", alpha = 0.7) +
  labs(
    title = "Distribution of % CG Methylation",
    x = "% Methylation per base",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# Classify CG windows based on methylation levels
unmethylated <- filter(cgsample, CG_methylation < 33)
hemimethylated <- filter(cgsample, CG_methylation >= 33 & CG_methylation <= 66)
methylated <- filter(cgsample, CG_methylation > 66)

# Convert categories to GRanges
unmethylated_GR <- GRanges(
  seqnames = unmethylated$chr, 
  ranges = IRanges(start = unmethylated$start, end = unmethylated$end), 
  CG_methylation = unmethylated$CG_methylation
)
hemimethylated_GR <- GRanges(
  seqnames = hemimethylated$chr, 
  ranges = IRanges(start = hemimethylated$start, end = hemimethylated$end), 
  CG_methylation = hemimethylated$CG_methylation
)
methylated_GR <- GRanges(
  seqnames = methylated$chr, 
  ranges = IRanges(start = methylated$start, end = methylated$end),
  CG_methylation = methylated$CG_methylation
)

# Plot histogram of CA methylation
ggplot(non_cg_windows, aes(x = CA_methylation)) +
  geom_histogram(binwidth = 5, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of % CA Methylation", 
    x = "% Methylation per base", 
    y = "Frequency"
  ) +
  theme_classic()

# Identify overlaps between CA windows and unmethylated CG regions
UM_over_hits <- findOverlaps(non_cg_windows_GR, unmethylated_GR)
UM <- cbind(
  unmethylated[subjectHits(UM_over_hits), ], 
  non_cg_windows[queryHits(UM_over_hits), ]
)
UMwindows <- UM[, c(1, 2, 3, 8, 9, 17, 18)]
UMwindows_GR_c <- makeGRangesFromDataFrame(
  UMwindows, 
  keep.extra.columns = TRUE, 
  start.field = "start", 
  end.field = "end",
  seqnames.field = "chr"
)

# Overlaps with hemimethylated CG windows
HM_over_hits <- findOverlaps(non_cg_windows_GR, hemimethylated_GR)
HM <- cbind(
  hemimethylated[subjectHits(HM_over_hits), ],
  non_cg_windows[queryHits(HM_over_hits), ]
)
HMwindows <- HM[, c(1, 2, 3, 8, 9, 17, 18)]
HMwindows_GR_c <- makeGRangesFromDataFrame(
  HMwindows, 
  keep.extra.columns = TRUE, 
  start.field = "start", 
  end.field = "end",
  seqnames.field = "chr"
)

# Overlaps with methylated CG windows
M_over_hits <- findOverlaps(non_cg_windows_GR, methylated_GR)
M <- cbind(
  methylated[subjectHits(M_over_hits), ],
  non_cg_windows[queryHits(M_over_hits), ]
)
Mwindows <- M[, c(1, 2, 3, 8, 9, 17, 18)]
Mwindows_GR_c <- makeGRangesFromDataFrame(
  Mwindows,
  keep.extra.columns = TRUE, 
  start.field = "start", 
  end.field = "end", 
  seqnames.field = "chr"
)

# Calculate fold change of CA methylation between methylated and unmethylated CG contexts
foldchange <- mean(Mwindows$CA_methylation) / mean(UMwindows$CA_methylation)

# Merge all categories into one data frame
cowindows <- rbind(
  UMwindows, 
  HMwindows, 
  Mwindows
  )

# Spearman correlation between CA and CG methylation
spear <- cor.test(
  cowindows$CA_methylation, 
  cowindows$CG_methylation, 
  method = "spearman"
)

# Linear regression using M-values
reg_mvalue <- summary(
  lm(mvalues_ca ~ mvalues_cg, data = cowindows)
)

# Summary statistics
resume <- data.frame(
  cg = mean(cgsample$CG_methylation),
  ca = mean(non_cg_windows_GR$CA_methylation),
  mca = mean(Mwindows$CA_methylation), 
  delta = mean(Mwindows$CA_methylation) - mean(non_cg_windows_GR$CA_methylation),
  FC = foldchange,
  Spearman = spear$estimate,
  slope = reg_mvalue$coefficients[2],
  R.adj = reg_mvalue$adj.r.squared,
  intercept = reg_mvalue$coefficients[1]
)

