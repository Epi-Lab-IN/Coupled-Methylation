##############################################################
# Gene-level aggregation and GO enrichment (g:Profiler)
# for CoM positive and negative regions
##############################################################

# -----------------------------
# Load required libraries
# -----------------------------
library(readr)
library(biomaRt)
library(GenomicRanges)
library(genomation)
library(dplyr)
library(gprofiler2)
library(ggplot2)
library(rrvgo)
library(org.Hs.eg.db)

# -----------------------------
# Input paths
# -----------------------------
path_COM_pos <- "/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/genbody_COM_pos_close10k.txt"
path_COM_neg <- "/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/genbody_COM_neg_close10k.txt"

# -----------------------------
# Load significant regions
# -----------------------------
sig_COM_pos <- read_delim(path_COM_pos, delim = "\t", show_col_types = FALSE)
sig_COM_neg <- read_delim(path_COM_neg, delim = "\t", show_col_types = FALSE)

# -----------------------------
# Get human gene annotation (protein-coding only)
# -----------------------------
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

all_genes <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl
)

# Clean chromosome names and convert to GRanges
all_genes <- all_genes %>%
  filter(grepl("^[0-9XY]+$", chromosome_name)) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name))

genes_gr <- makeGRangesFromDataFrame(
  all_genes,
  seqnames.field = "chromosome_name",
  start.field = "start_position",
  end.field = "end_position",
  keep.extra.columns = TRUE
)

# -----------------------------
# Helper function: Annotate and aggregate by gene
# -----------------------------
aggregate_by_gene <- function(sig_regions, genes_gr, diff_mode = c("max", "min")) {
  diff_mode <- match.arg(diff_mode)
  if (nrow(sig_regions) == 0) return(data.frame())
  
  sig_gr <- makeGRangesFromDataFrame(sig_regions, keep.extra.columns = TRUE)
  hits <- findOverlaps(sig_gr, genes_gr)
  if (length(hits) == 0) return(data.frame())
  
  df <- data.frame(
    ensembl_gene_id = genes_gr$ensembl_gene_id[subjectHits(hits)],
    gene_start      = start(genes_gr[subjectHits(hits)]),
    gene_end        = end(genes_gr[subjectHits(hits)]),
    diff            = sig_gr$diff[queryHits(hits)]
  )
  df$gene_size <- df$gene_end - df$gene_start + 1
  
  summarised <- df %>%
    group_by(ensembl_gene_id) %>%
    summarise(
      diff_stat    = if (diff_mode == "max") max(diff) else min(diff),
      n_regions    = n(),
      slope_per_kb = n_regions / (dplyr::first(gene_size) / 1000)
    ) %>%
    ungroup() %>%
    arrange(if (diff_mode == "max") dplyr::desc(diff_stat) else diff_stat)
  
  return(summarised)
}

# -----------------------------
# Apply to CoM positive and negative sets
# -----------------------------
genes_COM_pos <- aggregate_by_gene(sig_COM_pos, genes_gr, "max")
genes_COM_neg <- aggregate_by_gene(sig_COM_neg, genes_gr, "min")

# Filter genes with more than one significant region
genes_COM_pos_filt <- genes_COM_pos %>% filter(n_regions > 1)
genes_COM_neg_filt <- genes_COM_neg %>% filter(n_regions > 1)

# -----------------------------
# Prepare top genes for enrichment
# -----------------------------
top_genes_pos <- head(genes_COM_pos_filt$ensembl_gene_id, 1000)
top_genes_neg <- head(genes_COM_neg_filt$ensembl_gene_id, 1000)

# -----------------------------
# Functional enrichment (GO:BP)
# -----------------------------
gost_COM_pos <- gost(
  query = top_genes_pos,
  organism = "hsapiens",
  correction_method = "gSCS",
  sources = c("GO:BP")
)

gost_COM_neg <- gost(
  query = top_genes_neg,
  organism = "hsapiens",
  correction_method = "gSCS",
  sources = c("GO:BP")
)

# -----------------------------
# Helper function: Reduce redundant GO terms using semantic similarity
# -----------------------------
reduce_go_terms <- function(gost_result, group_name, top_n, threshold = 0.7, org_db = "org.Hs.eg.db") {
  if (is.null(gost_result$result)) return(data.frame())
  
  res <- gost_result$result %>%
    filter(term_size <= 500 & p_value < 0.05) %>%
    arrange(p_value)
  
  if (nrow(res) == 0) return(data.frame())
  
  # Calculate semantic similarity matrix
  sim_matrix <- calculateSimMatrix(
    res$term_id,
    orgdb = org_db,
    ont = "BP",
    method = "Rel"
  )
  
  # Use -log10(p_value) as score to prioritize most significant terms
  scores <- setNames(-log10(res$p_value), res$term_id)
  
  # Keep only scores for terms present in sim_matrix
  scores <- scores[names(scores) %in% rownames(sim_matrix)]
  
  # Reduce redundant terms
  reduced <- reduceSimMatrix(sim_matrix, scores, threshold = threshold, orgdb = org_db)
  
  # Keep only representative (parent) terms
  res_reduced <- res %>%
    filter(term_id %in% reduced$parent) %>%
    distinct(term_id, .keep_all = TRUE) %>%
    slice_head(n = top_n) %>%
    mutate(
      term_name  = factor(term_name, levels = rev(term_name)),
      gene_ratio = intersection_size / query_size,
      Group      = group_name
    )
  
  return(res_reduced)
}

# -----------------------------
# Extract and format top non-redundant terms
# -----------------------------
top_terms_pos <- reduce_go_terms(gost_COM_pos, "CoM +", top_n = 5, threshold = 0.7)
top_terms_neg <- reduce_go_terms(gost_COM_neg, "CoM -", top_n = 5, threshold = 0.7)

# Get unique selected terms across both groups
selected_terms <- unique(c(as.character(top_terms_pos$term_name),
                           as.character(top_terms_neg$term_name)))

# -----------------------------
# Expand results for selected terms (for plotting both groups together)
# -----------------------------
expand_terms <- function(gost_result, group_name, selected_terms) {
  if (is.null(gost_result$result)) return(data.frame())
  
  gost_result$result %>%
    filter(term_name %in% selected_terms) %>%
    mutate(
      gene_ratio = intersection_size / query_size,
      Group      = group_name
    )
}

plot_df <- bind_rows(
  expand_terms(gost_COM_pos, "CoM +", selected_terms),
  expand_terms(gost_COM_neg, "CoM -", selected_terms)
)

plot_df$Group     <- factor(plot_df$Group, levels = c("CoM +", "CoM -"))
plot_df$term_name <- factor(plot_df$term_name, levels = rev(unique(selected_terms)))

# -----------------------------
# Plot GO terms
# -----------------------------
ggplot(plot_df, aes(y = term_name, x = Group)) +
  geom_point(aes(size = gene_ratio, fill = p_value),
             shape = 21, color = "black", alpha = 0.8, stroke = 0.5) +
  scale_size(range = c(3, 8), name = "GeneRatio") +
  scale_fill_gradientn(
    colours = c("#7a0177", "#fbb4b9"),
    name = "p.adjust",
    trans = "log10",
    labels = scales::label_scientific()
  ) +
  theme_linedraw(base_size = 14) +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_text(size = 12),
    axis.ticks.x  = element_line(size = 1),
    axis.text.y   = element_text(size = 12),
    axis.ticks.y  = element_line(size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(size = 1, colour = "black"),
    legend.text   = element_text(size = 12),
    legend.title  = element_text(size = 13)
  ) +
  labs(y = NULL, x = NULL)

# -----------------------------
# Convert Ensembl IDs to gene symbols
# -----------------------------
annot_conv <- gconvert(
  query = c(genes_COM_pos_filt$ensembl_gene_id, genes_COM_neg_filt$ensembl_gene_id),
  organism = "hsapiens",
  target = "ENTREZGENE_ACC"
)

annot_conv <- annot_conv %>%
  dplyr::select(input, name, description) %>%
  dplyr::rename(ensembl_gene_id = input, gene_symbol = name, gene_description = description)

# Merge with aggregated results
genes_COM_pos_filt <- genes_COM_pos_filt %>%
  left_join(annot_conv, by = "ensembl_gene_id")

genes_COM_neg_filt <- genes_COM_neg_filt %>%
  left_join(annot_conv, by = "ensembl_gene_id")

# -----------------------------
# Save outputs
# -----------------------------
write.table(genes_COM_pos_filt, "genes_COM_pos_aggregated.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes_COM_neg_filt, "genes_COM_neg_aggregated.txt", sep = "\t", quote = FALSE, row.names = FALSE)