# ================================
# Libraries
# ================================
library(GenomicRanges)
library(genomation)
library(dplyr)
library(ggplot2)
library(readr)


# ================================
# Read data
# ================================
com <- read_delim("Escritorio/Animals/CHH/chimpance/Corr_with_age/CoM_neg_chimpance",
                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# Diccionario basado en mPanTro3-v2.1
map_refseq_to_chr <- c(
  "NC_072398.2" = "chr1",
  "NC_086015.1" = "chr2A",  
  "NC_072401.2" = "chr3",
  "NC_072402.2" = "chr4",
  "NC_072403.2" = "chr5",
  "NC_072404.2" = "chr6",
  "NC_072405.2" = "chr7",
  "NC_072406.2" = "chr8",
  "NC_072407.2" = "chr9",
  "NC_072408.2" = "chr10",
  "NC_072409.2" = "chr11",
  "NC_072410.2" = "chr12",
  "NC_072411.2" = "chr13",
  "NC_072412.2" = "chr14",
  "NC_072413.2" = "chr15",
  "NC_072414.2" = "chr16",
  "NC_072415.2" = "chr17",
  "NC_072416.2" = "chr18",
  "NC_072417.2" = "chr19",
  "NC_072418.2" = "chr20",
  "NC_072419.2" = "chr21",
  "NC_072420.2" = "chr22",
  "NC_086016.1" = "chr2B",
  "NC_072421.2" = "chrX",
  "NC_072422.2" = "chrY",
  "NC_001643.1" = "chrM"
)
com$chr <- map_refseq_to_chr[com$chr]
chromsizes <- read.table("/home/epi-lab/Escritorio/Database/panTro3.chrom.sizes.txt", header=FALSE) # We obtain from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
colnames(chromsizes) <- c("chr", "size")
chromsizes <- chromsizes[grepl("^chr[0-9XY]+$", chromsizes$chr), ]

# ================================
# Generate random windows proportional to chromosome sizes
# ================================
set.seed(42)
n_random <- nrow(com)
random_chrs <- sample(chromsizes$chr, n_random, replace = TRUE, prob = chromsizes$size / sum(chromsizes$size))
chrom_start_list <- setNames(lapply(chromsizes$size, function(sz) seq(1, sz - 100, by = 100)), chromsizes$chr)
random_starts <- mapply(function(chr, n) sample(chrom_start_list[[chr]], n, replace = TRUE),
                        chr = random_chrs, n = 1)
random_gr <- GRanges(seqnames = random_chrs, ranges = IRanges(start = random_starts, width = 100))
non_sig_sampled <- as.data.frame(random_gr)
colnames(non_sig_sampled)[1] <- "chr"


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

sig_df <- prepare_df(com, add_group = "CoM")
non_sig_df <- prepare_df(non_sig_sampled, add_group = "No correlation")

combined_df <- bind_rows(sig_df, non_sig_df)

# ================================
# Plot histogram comparing distances
# ================================
ggplot(combined_df, aes(x = distance_to_next, fill = group)) +
  geom_histogram(binwidth = 100, position = "identity", alpha = 1) +
  scale_fill_manual(values = c("CoM" = "#f05a3c", "No correlation" = "#babfba")) +
  labs(x = "Distance (bp)", y = "Frequency", fill = NULL) +
  theme_linedraw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 1.5, colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.ticks = element_line(size = 1),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  xlim(0, 20000) 


# ================================
# Función para obtener ventanas cercanas (<10kb) sin repeticiones
# ================================
get_close_windows <- function(df, max_dist = 10000) {
  gr <- GRanges(df$chr, IRanges(df$start, df$start + 99))
  
  # Encontramos todas las ventanas que están a menos de max_dist
  hits <- findOverlaps(gr, gr, maxgap = max_dist, type = "any")
  
  # Eliminamos auto-overlaps
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  
  # Tomamos todos los índices únicos de ventanas cercanas
  close_indices <- unique(c(queryHits(hits), subjectHits(hits)))
  
  # Devolvemos el dataframe filtrado sin repeticiones
  df[close_indices, ]
}
# ================================
# Aplicamos la función a los dos grupos
# ================================
close_sig_windows <- get_close_windows(sig_df, max_dist = 10002)
write.table(close_sig_windows[5:7], "com_pos_chimpance_clusters", sep = "\t", row.names = FALSE, col.names = FALSE)

close_non_sig_windows <- get_close_windows(non_sig_df, max_dist = 10002)

# ================================
# Fisher test: proporción de ventanas cercanas (<10kb)
# ================================

# Añadimos ID único a cada ventana
sig_df <- sig_df %>% mutate(id = paste0("sig_", row_number()))
non_sig_df <- non_sig_df %>% mutate(id = paste0("non_", row_number()))
df_all <- bind_rows(sig_df, non_sig_df)

# Detectamos las ventanas cercanas en cada grupo
close_sig <- get_close_windows(sig_df, max_dist = 10000)
close_non_sig <- get_close_windows(non_sig_df, max_dist = 10000)
random_gr <- makeGRangesFromDataFrame(close_non_sig_windows, keep.extra.columns = TRUE)
# Marcamos las que están en clusters cercanos
df_all <- df_all %>%
  mutate(close = id %in% c(close_sig$id, close_non_sig$id))

# Creamos tabla 2x2: grupo vs cercanía
contingency <- table(df_all$group, df_all$close)
print(contingency)

# Ejecutamos test de Fisher
fisher_result_exact <- fisher.test(contingency[2:1, ], alternative = "greater")
print(fisher_result_exact$p.value, digits=20)


# Optionally save results
write.table(close_sig_windows, "COM_pos_close10k_chimpance.txt", sep = "\t", row.names = FALSE)
com_gr <- makeGRangesFromDataFrame(close_sig_windows, keep.extra.columns = TRUE)

library(biomaRt)
library(GenomicRanges)

# 1. Conectar al mart de Chimpancé
chimp_ensembl <- useMart("ensembl", dataset = "ptroglodytes_gene_ensembl")

annot_chimp <- getBM(
  attributes = c(
    'ensembl_gene_id', 
    'chromosome_name', 
    'start_position', 
    'end_position', 
    'external_gene_name', 
    'hsapiens_homolog_ensembl_gene', 
    'hsapiens_homolog_orthology_type'
  ),
  mart = chimp_ensembl
)

# Filtrar para obtener solo los ortólogos más fiables (1 a 1)
orthos_clean <- annot_chimp[annot_chimp$hsapiens_homolog_orthology_type == "ortholog_one2one", ]

# Renombrar para mayor claridad
colnames(orthos_clean) <- c("gene_id_chimp", "chr", "start", "end", "gene_name", "gene_id_human", "orthology_type")
orthos_clean$chr <- paste0("chr", orthos_clean$chr)
genes_chimp_gr <- makeGRangesFromDataFrame(orthos_clean, keep.extra.columns = TRUE)


hits <- findOverlaps(com_gr, genes_chimp_gr)


genes<- cbind(orthos_clean[subjectHits(hits),], close_sig_windows[queryHits(hits),])
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

df_por_gen <- as.data.frame(df_por_gen)
genes_com_chimp <- df_por_gen[df_por_gen$n_ventanas>1,]


# 1. Solapamiento de las ventanas RANDOM con los genes de chimpancé
hits_random <- findOverlaps(random_gr, genes_chimp_gr)

# 2. Crear el dataframe de genes obtenidos por azar
genes_random_data <- cbind(orthos_clean[subjectHits(hits_random),], 
                           close_non_sig_windows[queryHits(hits_random),])
genes_random_data <- genes_random_data[-c(2:4)]
# 3. Resumir por gen (igual que hiciste con CoM)
genes_random_summary <- genes_random_data %>%
  group_by(gene_id_chimp) %>%
  summarise(
    n_ventanas = n(),
    gene_name = dplyr::first(gene_name),
    human_gene_id = dplyr::first(gene_id_human)
  ) %>%
  filter(n_ventanas > 1) # Mantener el mismo criterio de filtrado

#####
genes_human <-  read_delim("Escritorio/Resultados/CHH/genes/genes_COM_neg_aggregated.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

# IDs de genes humanos de referencia (tu agingMap)
human_ref_ids <- unique(genes_human$gene_symbol)

# Genes que coinciden en el grupo REAL
coincidencias_real <- sum(genes_com_chimp$gene_name %in% human_ref_ids)
total_real <- nrow(genes_com_chimp)

# Genes que coinciden en el grupo RANDOM
coincidencias_random <- sum(genes_random_summary$gene_name %in% human_ref_ids)
total_random <- nrow(genes_random_summary)


cont_matrix <- matrix(c(
  coincidencias_real, total_real - coincidencias_real,
  coincidencias_random, total_random - coincidencias_random
), nrow = 2, byrow = TRUE)

fisher_enrichment <- fisher.test(cont_matrix)

print(paste("P-value de enriquecimiento de ortólogos:", fisher_enrichment$p.value))
