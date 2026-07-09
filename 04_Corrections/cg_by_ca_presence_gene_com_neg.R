# ==============================================================================
# CG methylation analysis: CoM- genes vs random background genes,
# stratified by presence/absence of age-correlated (CA) regions
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
  com_file    = "Escritorio/Resultados/CHH/genes/genes_COM_neg_aggregated.txt",
  ca_file     = "Escritorio/Resultados/CHH/Corr_age/CA_bvalues_all.txt",
  cg_file     = "Escritorio/Resultados/CHH/Corr_age/CG_bvalues_all.txt",
  sample_cols = 4:10,     # columns holding per-sample beta values
  seed        = 123,
  bin_breaks  = c(-Inf, 33, 66, Inf),
  bin_labels  = c("Low\n(<33%)", "Mid\n(33-66%)", "High\n(>66%)")
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

classify_windows_cg <- function(cg_gr_full, genes_target_GR, ca_gr_ref,
                                group_label, bin_breaks, bin_labels) {
  
  # 1. CG windows inside the target genes (deduplicated)
  hits_genes  <- findOverlaps(cg_gr_full, genes_target_GR)
  cg_en_genes <- unique(cg_gr_full[queryHits(hits_genes)])
  
  # 2. Presence / absence of a CA region
  hits_ca    <- findOverlaps(cg_en_genes, ca_gr_ref)
  idx_con_ca <- unique(queryHits(hits_ca))
  
  b_vals  <- mcols(cg_en_genes)$b_value
  ca_flag <- ifelse(seq_along(b_vals) %in% idx_con_ca, "With CA", "Without CA")
  
  # 3. Fixed cutoffs at 33% / 66% (meaningful here because b_value is already a %)
  met_bin <- cut(b_vals, breaks = bin_breaks, labels = bin_labels, right = TRUE)
  
  message("-- ", group_label, " --")
  message("  N windows:     ", length(b_vals))
  message("  N with CA:     ", sum(ca_flag == "With CA"))
  message("  N without CA:  ", sum(ca_flag == "Without CA"))
  print(table(met_bin))
  
  data.frame(
    met       = b_vals,
    ca_status = ca_flag,
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

genes_GR <- harmonize_seqlevels(genes_GR, cg_gr)
ca_gr    <- harmonize_seqlevels(ca_gr, cg_gr)

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
random_GR <- harmonize_seqlevels(random_GR, cg_gr)

## ---- 5. Classify CG windows for both gene sets (single pass) ------------------


df_com  <- classify_windows_cg(cg_gr, genes_GR,  ca_gr, "CoM-",
                               CONFIG$bin_breaks, CONFIG$bin_labels)
df_rand <- classify_windows_cg(cg_gr, random_GR, ca_gr, "Random",
                               CONFIG$bin_breaks, CONFIG$bin_labels)

df_all <- rbind(df_com, df_rand)
df_all$group     <- factor(df_all$group,     levels = c("CoM-", "Random"))
df_all$ca_status <- factor(df_all$ca_status, levels = c("With CA", "Without CA"))
df_all$met_bin   <- factor(df_all$met_bin,   levels = CONFIG$bin_labels)

## ---- 6. Overall With-CA vs Without-CA comparison (per gene set) ---------------

overall_summary <- df_all %>%
  group_by(group, ca_status) %>%
  summarise(n = n(), mean_met = round(mean(met), 2),
            median_met = round(median(met), 2), .groups = "drop")
print(overall_summary)

overall_tests <- df_all %>%
  group_by(group) %>%
  summarise(
    mean_with    = round(mean(met[ca_status == "With CA"]),    2),
    mean_without = round(mean(met[ca_status == "Without CA"]), 2),
    p_wilcox = {
      x <- met[ca_status == "With CA"]
      y <- met[ca_status == "Without CA"]
      if (length(x) > 200000) x <- sample(x, 200000)
      if (length(y) > 200000) y <- sample(y, 200000)
      wilcox.test(x, y)$p.value
    },
    .groups = "drop"
  )
print(overall_tests)

## ---- 7. Boxplot: overall %mCG by CA status, faceted by gene set --------------

ggplot(df_all, aes(x = ca_status, y = met, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(50, 100)) +
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
  labs(y = "%mCG", x = NULL)

## ---- 8. Binned analysis: %mCG by methylation tercile x CA status -------------

binned_summary <- df_all %>%
  group_by(group, ca_status, met_bin) %>%
  summarise(n = n(), median = round(median(met), 2),
            mean = round(mean(met), 2), sd = round(sd(met), 2),
            .groups = "drop")
print(binned_summary)

pal <- c("With CA" = "#a961e8", "Without CA" = "#e68c4c")

ggplot(df_all, aes(x = met_bin, y = met, fill = ca_status)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.65) +
  coord_cartesian(ylim = c(0, 100)) +
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
  labs(y = "%mCG", x = NULL)

## ---- 9. Statistical tests within each methylation bin -------------------------
# Compares With-CA vs Without-CA methylation, separately for each gene set
# and each methylation tercile.

binned_tests_within_group <- df_all %>%
  group_by(group, met_bin) %>%
  summarise(
    p_wilcox = {
      x <- met[ca_status == "With CA"]
      y <- met[ca_status == "Without CA"]
      if (length(x) > 200000) x <- sample(x, 200000)
      if (length(y) > 200000) y <- sample(y, 200000)
      wilcox.test(x, y)$p.value
    },
    .groups = "drop"
  )
print(binned_tests_within_group)
