# ================================
# Libraries
# ================================
library(GenomicRanges)
library(genomation)
library(dplyr)
library(ggplot2)
library(readr)
library(biomaRt)

# ================================
# Configuration
# ================================
region_type <- "genbody"  # "genbody", "promoter", or "intergenic"

# Annotation filter function
region_filter <- function(class_matrix, type) {
  if(type == "genbody") {
    (class_matrix[,2] == 1 | class_matrix[,3] == 1) & class_matrix[,1] == 0
  } else if(type == "promoter") {
    class_matrix[,1] == 1
  } else if(type == "intergenic") {
    class_matrix[,2] == 0 & class_matrix[,3] == 0 & class_matrix[,1] == 0
  } else {
    stop("Unknown region type")
  }
}

# ================================
# Read data
# ================================
com <- read_delim("Escritorio/Animals/CHH/chimpance/Corr_with_age/CoM_pos_chimpance",
                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Diccionario RefSeq → chr
map_refseq_to_chr <- c(
  "NC_072398.2" = "chr1",  "NC_086015.1" = "chr2A", "NC_072401.2" = "chr3",
  "NC_072402.2" = "chr4",  "NC_072403.2" = "chr5",  "NC_072404.2" = "chr6",
  "NC_072405.2" = "chr7",  "NC_072406.2" = "chr8",  "NC_072407.2" = "chr9",
  "NC_072408.2" = "chr10", "NC_072409.2" = "chr11", "NC_072410.2" = "chr12",
  "NC_072411.2" = "chr13", "NC_072412.2" = "chr14", "NC_072413.2" = "chr15",
  "NC_072414.2" = "chr16", "NC_072415.2" = "chr17", "NC_072416.2" = "chr18",
  "NC_072417.2" = "chr19", "NC_072418.2" = "chr20", "NC_072419.2" = "chr21",
  "NC_072420.2" = "chr22", "NC_086016.1" = "chr2B", "NC_072421.2" = "chrX",
  "NC_072422.2" = "chrY",  "NC_001643.1" = "chrM"
)
com$chr <- map_refseq_to_chr[com$chr]

chromsizes <- read.table("/home/epi-lab/Escritorio/Database/panTro3.chrom.sizes.txt", header = FALSE)
colnames(chromsizes) <- c("chr", "size")
chromsizes <- chromsizes[grepl("^chr[0-9XY]+$", chromsizes$chr), ]

# ================================
# Annotate and filter by genomic region
# ================================
gene_features <- readTranscriptFeatures("/home/epi-lab/Escritorio/Animals/chimpance/hgTables.bed",
                                        remove.unusual = TRUE)

annotate_and_filter <- function(df, region_type) {
  gr         <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  annotation <- annotateWithGeneParts(gr, gene_features)
  class_mat  <- getMembers(annotation)
  df[region_filter(class_mat, region_type), ]
}

com_filtered <- annotate_and_filter(com, region_type)
cat("CoM windows after region filter:", nrow(com_filtered), "\n")

# ================================
# Generate random windows proportional to chromosome sizes
# ================================

CA <- read_delim("~/Escritorio/Animals/CHH/chimpance/chimpance_ca_comunes_bvalues.tsv", delim = "\t", trim_ws = TRUE)

CG <- read_delim("~/Escritorio/Animals/CHH/chimpance/chimpance_cg_comunes_bvalues.tsv", delim = "\t", trim_ws = TRUE)
comun <- rbind(CA[1:3], CG[1:3]) %>% distinct()
comun$chr <- map_refseq_to_chr[comun$chr]
# Samplear del universo filtrado
set.seed(42)
n_random <- nrow(com)
random_chrs <- sample(chromsizes$chr, n_random, replace = TRUE, prob = chromsizes$size / sum(chromsizes$size))
chrom_start_list <- setNames(lapply(chromsizes$size, function(sz) seq(1, sz - 100, by = 100)), chromsizes$chr)
random_starts <- mapply(function(chr, n) sample(chrom_start_list[[chr]], n, replace = TRUE),
                        chr = random_chrs, n = 1)
random_gr <- GRanges(seqnames = random_chrs, ranges = IRanges(start = random_starts, width = 100))
non_sig_sampled <- as.data.frame(random_gr)
colnames(non_sig_sampled)[1] <- "chr"


non_sig_filtered <- annotate_and_filter(non_sig_sampled, region_type)

# ================================
# Prepare dataframes
# ================================
prepare_df <- function(df, add_group = NULL) {
  df %>%
    arrange(chr, start) %>%
    group_by(chr) %>%
    mutate(distance_to_next = lead(start) - end) %>%
    ungroup() %>%
    {if(!is.null(add_group)) mutate(., group = add_group) else .}
}

sig_df     <- prepare_df(com_filtered,  add_group = "CoM")
non_sig_df <- prepare_df(non_sig_filtered, add_group = "No correlation")
combined_df <- bind_rows(sig_df, non_sig_df)

# ================================
# Plot histogram
# ================================
ggplot(combined_df, aes(x = distance_to_next, fill = group)) +
  geom_histogram(binwidth = 100, position = "identity", alpha = 1) +
  scale_fill_manual(values = c("CoM" = "#f05a3c", "No correlation" = "#babfba")) +
  labs(x = "Distance (bp)", y = "Frequency", fill = NULL) +
  theme_linedraw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(size = 1.5, colour = "black"),
    axis.text        = element_text(size = 18),
    axis.title       = element_text(size = 18),
    axis.ticks       = element_line(size = 1),
    legend.text      = element_text(size = 16),
    legend.title     = element_text(size = 18)
  ) +
  xlim(0, 20000)

# ================================
# Close windows (<10kb)
# ================================
get_close_windows <- function(df, max_dist = 10000) {
  gr   <- GRanges(df$chr, IRanges(df$start, df$start + 99))
  hits <- findOverlaps(gr, gr, maxgap = max_dist, type = "any")
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  close_indices <- unique(c(queryHits(hits), subjectHits(hits)))
  df[close_indices, ]
}

close_sig_windows     <- get_close_windows(sig_df,     max_dist = 10002)
close_non_sig_windows <- get_close_windows(non_sig_df, max_dist = 10002)

# ================================
# Fisher test
# ================================
sig_df     <- sig_df     %>% mutate(id = paste0("sig_",  row_number()))
non_sig_df <- non_sig_df %>% mutate(id = paste0("non_",  row_number()))
df_all     <- bind_rows(sig_df, non_sig_df)

close_sig     <- get_close_windows(sig_df,     max_dist = 10000)
close_non_sig <- get_close_windows(non_sig_df, max_dist = 10000)

random_gr <- makeGRangesFromDataFrame(close_non_sig_windows, keep.extra.columns = TRUE)

df_all <- df_all %>%
  mutate(close = id %in% c(close_sig$id, close_non_sig$id))

contingency  <- table(df_all$group, df_all$close)
print(contingency)

fisher_result_exact <- fisher.test(contingency[2:1, ], alternative = "greater")
print(fisher_result_exact$p.value, digits = 20)

# ================================
# Gene annotation (biomaRt)
# ================================
com_gr <- makeGRangesFromDataFrame(close_sig_windows, keep.extra.columns = TRUE)

chimp_ensembl <- useMart("ensembl", dataset = "ptroglodytes_gene_ensembl")

annot_chimp <- getBM(
  attributes = c(
    'ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position',
    'external_gene_name', 'hsapiens_homolog_ensembl_gene',
    'hsapiens_homolog_orthology_type'
  ),
  mart = chimp_ensembl
)

orthos_clean <- annot_chimp[annot_chimp$hsapiens_homolog_orthology_type == "ortholog_one2one", ]
colnames(orthos_clean) <- c("gene_id_chimp", "chr", "start", "end",
                            "gene_name", "gene_id_human", "orthology_type")
orthos_clean$chr      <- paste0("chr", orthos_clean$chr)
genes_chimp_gr        <- makeGRangesFromDataFrame(orthos_clean, keep.extra.columns = TRUE)

hits  <- findOverlaps(com_gr, genes_chimp_gr)
genes <- cbind(orthos_clean[subjectHits(hits), ], close_sig_windows[queryHits(hits), ])
genes <- genes[-c(2:4)]

df_por_gen <- genes %>%
  group_by(gene_id_chimp) %>%
  summarise(
    n_ventanas    = n(),
    R             = mean(R.ages),
    p.value       = mean(p.value),
    slope         = mean(slope.ages),
    gene_name     = dplyr::first(gene_name),
    gene_id_chimp = dplyr::first(gene_id_chimp),
    human_gene_id = dplyr::first(gene_id_human),
    chr           = dplyr::first(chr),
    start_gen     = min(start),
    end_gen       = max(end)
  ) %>%
  ungroup()

df_por_gen     <- as.data.frame(df_por_gen)
genes_com_chimp <- df_por_gen[df_por_gen$n_ventanas > 1, ]

# ================================
# Random vs genes
# ================================
hits_random        <- findOverlaps(random_gr, genes_chimp_gr)
genes_random_data  <- cbind(orthos_clean[subjectHits(hits_random), ],
                            close_non_sig_windows[queryHits(hits_random), ])
genes_random_data  <- genes_random_data[-c(2:4)]

genes_random_summary <- genes_random_data %>%
  group_by(gene_id_chimp) %>%
  summarise(
    n_ventanas    = n(),
    gene_name     = dplyr::first(gene_name),
    human_gene_id = dplyr::first(gene_id_human)
  ) %>%
  filter(n_ventanas > 1)

# ================================
# Enrichment Fisher test vs human aging map
# ================================
genes_human <- read_delim("Escritorio/Resultados/CHH/genes/genes_COM_pos_aggregated.txt",
                          delim = "\t", escape_double = FALSE, trim_ws = TRUE)

human_ref_ids       <- unique(genes_human$gene_symbol)
coincidencias_real  <- sum(genes_com_chimp$gene_name   %in% human_ref_ids)
total_real          <- nrow(genes_com_chimp)
coincidencias_random <- sum(genes_random_summary$gene_name %in% human_ref_ids)
total_random        <- nrow(genes_random_summary)

cont_matrix <- matrix(c(
  coincidencias_real,   total_real   - coincidencias_real,
  coincidencias_random, total_random - coincidencias_random
), nrow = 2, byrow = TRUE)

fisher_enrichment <- fisher.test(cont_matrix)
print(paste("P-value de enriquecimiento de ortólogos:", fisher_enrichment$p.value))






# Background = todos los genes con ortólogo humano
total_chimp_genes  <- nrow(orthos_clean)
coincidencias_background <- sum(orthos_clean$gene_name %in% human_ref_ids)

cont_matrix <- matrix(c(
  coincidencias_real,          total_real - coincidencias_real,
  coincidencias_background,    total_chimp_genes - coincidencias_background
), nrow = 2, byrow = TRUE)

fisher_enrichment <- fisher.test(cont_matrix)


cat("Hypergeometric enrichment test:\n")
cat("  CoM+ genes overlapping human aging genes:", coincidencias_real, "/", total_real, "\n")
cat("  Expected by chance:", round(total_real * length(human_ref_ids) / total_chimp_genes, 2), "\n")
cat("  p-value:", phyper(
  q = coincidencias_real - 1,
  m = length(human_ref_ids),
  n = total_chimp_genes - length(human_ref_ids),
  k = total_real,
  lower.tail = FALSE
), "\n")


n_perm <- 1000
overlaps_perm <- replicate(n_perm, {
  genes_random <- sample(orthos_clean$gene_name, total_real)
  sum(genes_random %in% human_ref_ids)
})

# Z-score como estadístico
z_score <- (coincidencias_real - mean(overlaps_perm)) / sd(overlaps_perm)

# P-value empírico
p_perm <- mean(overlaps_perm >= coincidencias_real)

cat("Permutation test:\n")
cat("  Observed overlaps:  ", coincidencias_real, "\n")
cat("  Expected (mean):    ", round(mean(overlaps_perm), 2), "\n")
cat("  SD:                 ", round(sd(overlaps_perm), 2), "\n")
cat("  Z-score:            ", round(z_score, 2), "\n")
cat("  p-value (empirical):", p_perm, "\n")

genes_comun <- genes_com_chimp[genes_com_chimp$gene_name   %in% human_ref_ids,]
write.table(genes_comun, file = "genes_com_pos_comun_chimp_human.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
