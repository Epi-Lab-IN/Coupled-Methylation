# ===============================================================
# Analysis of CA methylation in genes CoM+
# ===============================================================

# -------------------------
# 0. Load required libraries
# -------------------------
library(biomaRt)
library(GenomicRanges)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------
# 1. Load CoM+ regions
# -------------------------
# Load genomic regions close to CoM+ (positive slope) elements
com_pos_regions <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/genbody_COM_pos_close10k.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)
com_pos_regions$chr <- gsub("chr", "", com_pos_regions$chr)
# -------------------------
# 2. Retrieve Ensembl gene annotations
# -------------------------
# Connect to Ensembl BioMart (human dataset)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve coordinates of protein-coding genes
protein_coding_genes <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl
)

# Convert gene table to GRanges
genes_GR <- makeGRangesFromDataFrame(
  protein_coding_genes,
  keep.extra.columns = TRUE,
  start.field = "start_position",
  end.field = "end_position",
  seqnames.field = "chromosome_name"
)

# -------------------------
# 3. Find overlaps between CoM+ regions and genes
# -------------------------
# Convert input to GRanges and find overlaps
overlap_hits <- findOverlaps(as(com_pos_regions, "GRanges"), genes_GR)

# Combine overlapping genes with their corresponding CoM+ regions
overlapping_genes <- cbind(
  protein_coding_genes[subjectHits(overlap_hits), ],
  com_pos_regions[queryHits(overlap_hits), ]
)

# -------------------------
# 4. Load CA methylation windows
# -------------------------
ca_windows <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/CHH/EMseq_ca_comunes_bvalues.tsv",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Rename CA methylation columns (one per age)
colnames(ca_windows)[4:10] <- c(
  "y35_ca", "y45_ca", "y55_ca", "y65_ca",
  "y75_ca", "y85_ca", "y103_ca"
)

# Convert CA windows to GRanges
ca_windows_GR <- makeGRangesFromDataFrame(ca_windows, keep.extra.columns = TRUE)

# -------------------------
# 5. Define list of genes to visualize
# -------------------------
# Example gene (replace or extend this list)
target_genes <- c("ENSG00000179242", "ENSG00000183765", "ENSG00000146648", "ENSG00000141736")

# -------------------------
# 6. Function: process one gene and compute mean mCA per age
# -------------------------
process_gene <- function(gene_id, overlapping_genes, ca_windows_GR) {
  
  # Select all CoM+ regions linked to this gene
  gene_regions <- unique(overlapping_genes[overlapping_genes$ensembl_gene_id == gene_id, ])
  
  # Convert selected regions to GRanges
  gene_regions_GR <- makeGRangesFromDataFrame(
    gene_regions,
    keep.extra.columns = TRUE,
    start.field = "start",
    end.field = "end",
    seqnames.field = "chr"
  )
  
  # Find overlapping CA methylation windows
  hits <- findOverlaps(gene_regions_GR, ca_windows_GR)
  
  # If no overlaps, skip this gene
  if (length(hits) == 0) return(NULL)
  
  # Combine coordinates and methylation data
  df_ca <- cbind(
    as.data.frame(gene_regions_GR[queryHits(hits)]),
    as.data.frame(mcols(ca_windows_GR[subjectHits(hits)]))
  )
  
  # Keep only methylation columns (yXX_ca)
  df_meth <- df_ca %>%
    dplyr::select(starts_with("y"))
  
  # Collapse multiple windows → mean methylation per age
  df_gene <- df_meth %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    mutate(gene_id = gene_id)
  
  return(df_gene)
}

# -------------------------
# 7. Apply function to all selected genes
# -------------------------
genes_methylation_df <- lapply(
  target_genes,
  process_gene,
  overlapping_genes = overlapping_genes,
  ca_windows_GR = ca_windows_GR
) %>%
  bind_rows()

# -------------------------
# 8. Reshape data for heatmap visualization
# -------------------------
# Convert wide format (columns per age) to long format
df_long <- genes_methylation_df %>%
  pivot_longer(
    cols = starts_with("y"),
    names_to = "sample",
    values_to = "methyl_CA"
  )

# Create a readable label for each age
sample_labels <- c(
  "y35_ca"  = "35 years",
  "y45_ca"  = "45 years",
  "y55_ca"  = "55 years",
  "y65_ca"  = "65 years",
  "y75_ca"  = "75 years",
  "y85_ca"  = "85 years",
  "y103_ca" = "103 years"
)

# Apply factor levels to keep age order
df_long$sample <- factor(
  sample_labels[df_long$sample],
  levels = c("35 years", "45 years", "55 years", "65 years", "75 years", "85 years", "103 years")
)

# Ensure genes are plotted in the specified order
df_long$gene_id <- factor(df_long$gene_id, levels = rev(target_genes))

# -------------------------
# 9. Plot heatmap (one row per gene, one column per age)
# -------------------------
ggplot(df_long, aes(x = sample, y = gene_id, fill = methyl_CA)) +
  geom_tile(width = 1, height = 1, alpha = 0.8) +
  scale_fill_gradient(low = "white", high = "#0072B2", name = "% mCA") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
##################################################################################################
# ===============================================================
# Analysis of CA methylation in genes CoM+
# ===============================================================

# -------------------------
# 0. Load required libraries
# -------------------------
library(biomaRt)
library(GenomicRanges)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------
# 1. Load CoM- regions
# -------------------------
# Load genomic regions close to CoM+ (positive slope) elements
com_neg_regions <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/genbody_COM_neg_close10k.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

com_neg_regions$chr <- gsub("chr", "", com_neg_regions$chr)
# -------------------------
# 3. Find overlaps between CoM- regions and genes
# -------------------------
# Convert input to GRanges and find overlaps
overlap_hits <- findOverlaps(as(com_neg_regions, "GRanges"), genes_GR)

# Combine overlapping genes with their corresponding CoM+ regions
overlapping_genes <- cbind(
  protein_coding_genes[subjectHits(overlap_hits), ],
  com_neg_regions[queryHits(overlap_hits), ]
)

# -------------------------
# 4. Load CGmethylation windows
# -------------------------
cg_windows <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/EMseq_cg_comunes_bvalues.tsv",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Rename CG methylation columns (one per age)
colnames(cg_windows)[4:10] <- c(
  "y35_cg", "y45_cg", "y55_cg", "y65_cg",
  "y75_cg", "y85_cg", "y103_cg"
)

# Convert CA windows to GRanges
cg_windows_GR <- makeGRangesFromDataFrame(cg_windows, keep.extra.columns = TRUE)

# -------------------------
# 5. Define list of genes to visualize
# -------------------------
# Example gene (replace or extend this list)
target_genes <- c("ENSG00000142192", "ENSG00000048740", "ENSG00000119912", "ENSG00000116731")

# -------------------------
# 7. Apply function to all selected genes
# -------------------------
genes_methylation_df <- lapply(
  target_genes,
  process_gene,
  overlapping_genes = overlapping_genes,
  ca_windows_GR = cg_windows_GR
) %>%
  bind_rows()

# -------------------------
# 8. Reshape data for heatmap visualization
# -------------------------
# Convert wide format (columns per age) to long format
df_long <- genes_methylation_df %>%
  pivot_longer(
    cols = starts_with("y"),
    names_to = "sample",
    values_to = "methyl_CG"
  )

# Create a readable label for each age
sample_labels <- c(
  "y35_cg"  = "35 years",
  "y45_cg"  = "45 years",
  "y55_cg"  = "55 years",
  "y65_cg"  = "65 years",
  "y75_cg"  = "75 years",
  "y85_cg"  = "85 years",
  "y103_cg" = "103 years"
)

# Apply factor levels to keep age order
df_long$sample <- factor(
  sample_labels[df_long$sample],
  levels = c("35 years", "45 years", "55 years", "65 years", "75 years", "85 years", "103 years")
)

# Ensure genes are plotted in the specified order
df_long$gene_id <- factor(df_long$gene_id, levels = rev(target_genes))

# -------------------------
# 9. Plot heatmap (one row per gene, one column per age)
# -------------------------
ggplot(df_long, aes(x = sample, y = gene_id, fill = methyl_CG)) +
  geom_tile(width = 1, height = 1, alpha = 0.8) +
  scale_fill_gradient(low = "white", high = "#C22506", name = "% mCG") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
