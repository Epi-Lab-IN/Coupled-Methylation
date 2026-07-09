# ==============================================================================
# CA methylation analysis: CoM+ genes vs random background genes,
# stratified by presence/absence of CG regions
# ==============================================================================

## ---- 0. Libraries -----------------------------------------------------------
library(readr)
library(biomaRt)
library(GenomicRanges)
library(ggplot2)
library(dplyr)

## ---- 1. Configuration -------------------------------------------------------
# Centralize file paths and parameters here so they only need to be edited once
CONFIG <- list(
  com_file    = "Escritorio/Resultados/CHH/genes/genes_COM_pos_aggregated.txt",
  ca_file     = "Escritorio/Resultados/CHH/Corr_age/CA_bvalues_all.txt",
  cg_file     = "Escritorio/Resultados/CHH/Corr_age/CG_bvalues_all.txt",
  sample_cols = 4:10,     # columns holding per-sample beta values
  seed        = 123,
  bin_probs   = c(0.33, 0.66),  # quantiles used to split into Low/Mid/High
  bin_labels  = c("Low\n(<Q33)", "Mid\n(Q33-Q66)", "High\n(>Q66)")
)

## ---- 2. Helper functions -----------------------------------------------------

#' Read a beta-value table and add a mean b_value column
read_bval_table <- function(path, sample_cols) {
  df <- read_delim(path, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                   show_col_types = FALSE)
  df$b_value <- rowMeans(df[, sample_cols], na.rm = TRUE)
  df
}

#' Download protein-coding gene coordinates from Ensembl
get_protein_coding_genes <- function(dataset = "hsapiens_gene_ensembl") {
  ensembl <- useEnsembl(biomart = "genes", dataset = dataset)
  getBM(
    # Fixed: the original had a stray empty string "" in this vector, which is
    # not a valid attribute and would break/confuse the query.
    attributes = c('ensembl_gene_id', 'chromosome_name',
                   'start_position', 'end_position'),
    filters = 'biotype',
    values  = "protein_coding",
    mart    = ensembl
  )
}

harmonize_seqlevels <- function(gr_to_fix, gr_reference) {
  tryCatch({
    seqlevelsStyle(gr_to_fix) <- seqlevelsStyle(gr_reference)
  }, error = function(e) {
    message("Could not harmonize seqlevels style automatically: ", conditionMessage(e))
  })
  gr_to_fix
}


classify_windows <- function(ca_gr_full, genes_target_GR, cg_gr_ref,
                             group_label, bin_probs, bin_labels) {
  
  # 1. CA windows inside the target genes (deduplicated)
  hits_genes  <- findOverlaps(ca_gr_full, genes_target_GR)
  ca_en_genes <- unique(ca_gr_full[queryHits(hits_genes)])
  
  # 2. Presence / absence of a CG region
  hits_cg    <- findOverlaps(ca_en_genes, cg_gr_ref)
  idx_con_cg <- unique(queryHits(hits_cg))
  
  b_vals  <- mcols(ca_en_genes)$b_value
  cg_flag <- ifelse(seq_along(b_vals) %in% idx_con_cg, "With CG", "Without CG")
  
  # 3. Group-specific quantiles (computed over all CA windows of this group)
  q33 <- quantile(b_vals, bin_probs[1], na.rm = TRUE)
  q66 <- quantile(b_vals, bin_probs[2], na.rm = TRUE)
  
  message("-- ", group_label, " --")
  message("  Q33: ", round(q33, 2), "%")
  message("  Q66: ", round(q66, 2), "%")
  message("  N windows:    ", length(b_vals))
  message("  N with CG:    ", sum(cg_flag == "With CG"))
  message("  N without CG: ", sum(cg_flag == "Without CG"))
  
  # 4. Bin classification using the group's own quantiles
  met_bin <- cut(b_vals, breaks = c(-Inf, q33, q66, Inf),
                 labels = bin_labels, right = TRUE)
  
  data.frame(
    met       = b_vals,
    cg_status = cg_flag,
    met_bin   = met_bin,
    group     = group_label
  )
}

## ---- 3. Load data (only once) ------------------------------------------------


com <- read_delim(CONFIG$com_file, delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE, show_col_types = FALSE)

all_genes <- get_protein_coding_genes()

genes <- merge(com, all_genes, by = "ensembl_gene_id")

genes_GR <- makeGRangesFromDataFrame(
  genes, keep.extra.columns = TRUE,
  start.field = "start_position", end.field = "end_position",
  seqnames.field = "chromosome_name"
)

ca    <- read_bval_table(CONFIG$ca_file, CONFIG$sample_cols)
ca_gr <- makeGRangesFromDataFrame(ca, keep.extra.columns = TRUE)

cg    <- read_bval_table(CONFIG$cg_file, CONFIG$sample_cols)
cg_gr <- makeGRangesFromDataFrame(cg, keep.extra.columns = TRUE)


genes_GR <- harmonize_seqlevels(genes_GR, ca_gr)
cg_gr    <- harmonize_seqlevels(cg_gr, ca_gr)

## ---- 4. Random background gene set -------------------------------------------

n_genes_com      <- nrow(genes)
background_genes <- all_genes[!(all_genes$ensembl_gene_id %in% genes$ensembl_gene_id), ]

set.seed(CONFIG$seed)
random_genes <- background_genes[sample(nrow(background_genes), n_genes_com), ]

random_GR <- makeGRangesFromDataFrame(
  random_genes, keep.extra.columns = TRUE,
  start.field = "start_position", end.field = "end_position",
  seqnames.field = "chromosome_name"
)
random_GR <- harmonize_seqlevels(random_GR, ca_gr)

## ---- 5. Classify CA windows for both gene sets (single pass) ------------------


df_com  <- classify_windows(ca_gr, genes_GR,  cg_gr, "CoM+",
                            CONFIG$bin_probs, CONFIG$bin_labels)
df_rand <- classify_windows(ca_gr, random_GR, cg_gr, "Random",
                            CONFIG$bin_probs, CONFIG$bin_labels)

df_all <- rbind(df_com, df_rand)
df_all$group     <- factor(df_all$group,     levels = c("CoM+", "Random"))
df_all$cg_status <- factor(df_all$cg_status, levels = c("With CG", "Without CG"))
df_all$met_bin   <- factor(df_all$met_bin,   levels = CONFIG$bin_labels)

## ---- 6. Overall With-CG vs Without-CG comparison (per gene set) ---------------

overall_summary <- df_all %>%
  group_by(group, cg_status) %>%
  summarise(n = n(), mean_met = round(mean(met), 2),
            median_met = round(median(met), 2), .groups = "drop")
print(overall_summary)

overall_tests <- df_all %>%
  group_by(group) %>%
  summarise(
    mean_with    = round(mean(met[cg_status == "With CG"]),    2),
    mean_without = round(mean(met[cg_status == "Without CG"]), 2),
    p_wilcox = {
      x <- met[cg_status == "With CG"]
      y <- met[cg_status == "Without CG"]
      if (length(x) > 200000) x <- sample(x, 200000)
      if (length(y) > 200000) y <- sample(y, 200000)
      wilcox.test(x, y)$p.value
    },
    .groups = "drop"
  )
print(overall_tests)

## ---- 7. Boxplot: overall %mCA by CG status, faceted by gene set --------------

ggplot(df_all, aes(x = cg_status, y = met, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 30)) +
  facet_wrap(~ group) +
  theme_linedraw(base_size = 14) +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 12),
    axis.ticks.x     = element_line(linewidth = 1),
    axis.text.y      = element_text(size = 12),
    axis.ticks.y     = element_line(linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(linewidth = 1, colour = "black"),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 13)
  ) +
  labs(y = "%mCA", x = NULL)

## ---- 8. Binned analysis: %mCA by group-specific tercile x CG status ----------


binned_summary <- df_all %>%
  group_by(group, cg_status, met_bin) %>%
  summarise(n = n(), median = round(median(met), 2),
            mean = round(mean(met), 2), sd = round(sd(met), 2),
            .groups = "drop")
print(binned_summary)

pal <- c("With CG" = "#a961e8", "Without CG" = "#e68c4c")

ggplot(df_all, aes(x = met_bin, y = met, fill = cg_status)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.65) +
  coord_cartesian(ylim = c(0, 40)) +
  facet_wrap(~ group, ncol = 2) +
  scale_fill_manual(values = pal, name = "Gene set") +
  theme_linedraw(base_size = 14) +
  theme(
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 10),
    axis.text.y      = element_text(size = 12),
    axis.ticks       = element_line(linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(linewidth = 1, colour = "black"),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(colour = "white", size = 13, face = "bold"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 13)
  ) +
  labs(y = "%mCA", x = NULL)

## ---- 9. Statistical tests within each methylation bin -------------------------
# Compares With-CG vs Without-CG methylation, separately for each gene set
# and each methylation tercile.

binned_tests_within_group <- df_all %>%
  group_by(group, met_bin) %>%
  summarise(
    p_wilcox = {
      x <- met[cg_status == "With CG"]
      y <- met[cg_status == "Without CG"]
      if (length(x) > 200000) x <- sample(x, 200000)
      if (length(y) > 200000) y <- sample(y, 200000)
      wilcox.test(x, y)$p.value
    },
    .groups = "drop"
  )
print(binned_tests_within_group)
