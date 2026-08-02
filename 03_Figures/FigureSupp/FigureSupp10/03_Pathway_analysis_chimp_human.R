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

# -----------------------------
# Input paths
# -----------------------------
path_COM_pos <- "/home/epi-lab/genes_com_pos_chimp.tsv"
path_COM_neg <- "/home/epi-lab/Escritorio/Resultados/CHH/genes/genes_COM_pos_aggregated.txt"

# -----------------------------
# Load significant regions
# -----------------------------
sig_COM_pos <- read_delim(path_COM_pos, delim = "\t", show_col_types = FALSE)
sig_COM_neg <- read_delim(path_COM_neg, delim = "\t", show_col_types = FALSE)

# -----------------------------
# Prepare genes for enrichment
# -----------------------------

# -----------------------------
# Functional enrichment (GO:BP)
# -----------------------------
gost_COM_pos <- gost(
  query             = sig_COM_pos$human_gene_id,
  organism          = "hsapiens",
  correction_method = "gSCS",
  sources           = "GO:BP"
)

gost_COM_neg <- gost(
  query             = sig_COM_neg$ensembl_gene_id,
  organism          = "hsapiens",
  correction_method = "gSCS",
  sources           = "GO:BP"
)

# -----------------------------
# Extract and format top terms
# -----------------------------
get_top_terms <- function(gost_result, group_name, top_n = 5) {
  if (is.null(gost_result$result)) return(data.frame())
  
  gost_result$result %>%
    dplyr::arrange(p_value) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      gene_ratio = intersection_size / query_size,
      Group      = group_name
    )
}

top_terms_pos <- get_top_terms(gost_COM_pos, "Chimpanzee", 10)
top_terms_neg <- get_top_terms(gost_COM_neg, "Human", 10)

selected_terms <- unique(c(top_terms_pos$term_name, top_terms_neg$term_name))

# -----------------------------
# Expand all results for plotting
# -----------------------------
expand_terms <- function(gost_result, group_name, selected_terms) {
  if (is.null(gost_result$result)) return(data.frame())
  
  gost_result$result %>%
   # filter(term_size <= 5000 & p_value < 0.05) %>%
    dplyr::filter(term_name %in% selected_terms) %>%
    dplyr::mutate(
      gene_ratio = intersection_size / query_size,
      Group      = group_name
    )
}

plot_df <- bind_rows(
  expand_terms(gost_COM_pos, "Chimpanzee", selected_terms),
  expand_terms(gost_COM_neg, "Human", selected_terms)
)

plot_df$Group     <- factor(plot_df$Group,     levels = c("Chimpanzee", "Human"))
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
    name    = "p.adjust",
    trans   = "log10",
    labels  = scales::label_scientific()
  ) +
  theme_linedraw(base_size = 14) +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 12),
    axis.ticks.x     = element_line(size = 1),
    axis.text.y      = element_text(size = 12),
    axis.ticks.y     = element_line(size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(size = 1, colour = "black"),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 13)
  ) +
  labs(y = NULL, x = NULL)

# -----------------------------
# Save outputs
# -----------------------------
write.table(genes_COM_pos_filt, "genes_COM_pos_aggregated.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes_COM_neg_filt, "genes_COM_neg_aggregated.txt", sep = "\t", quote = FALSE, row.names = FALSE)