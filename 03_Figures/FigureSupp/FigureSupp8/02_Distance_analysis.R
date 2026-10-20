# ================================
# Libraries
# ================================
library(GenomicRanges)
library(genomation)
library(dplyr)
library(ggplot2)
library(readr)

# ================================
# Parameters for analysis
# ================================
# Define the region type to analyze: "genbody", "promoter", "intergenic"
region_type <- "genbody"

# Annotation filter function based on region type
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
com <- read_delim("Escritorio/Resultados/Corr_with_age/COM_pos_filtered",
                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

chromsizes <- read.table("/home/epi-lab/Escritorio/Database/hg38.chrom.sizes", header=FALSE) # We obtain from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
colnames(chromsizes) <- c("chr", "size")
chromsizes <- chromsizes[grepl("^chr[0-9XY]+$", chromsizes$chr), ]

gene_features <- readTranscriptFeatures("/home/epi-lab/Escritorio/Database/hg38.bed", remove.unusual = TRUE) 
#We obtain bed file from https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=3856151035_MYca0IwU34JARJDrh0lGVmCABp39. Selecting hg38 genome, group: genes and gene predictions, Track: NCBI Refseq all

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
non_sig_sampled$chr <- gsub("^chr", "", non_sig_sampled$chr)

# ================================
# Annotate significant and non-significant windows
# ================================
annotate_and_filter <- function(df, gr, region_type) {
  annotation <- annotateWithGeneParts(as(gr, "GRanges"), gene_features)
  class_mat <- getMembers(annotation)
  df[region_filter(class_mat, region_type), ]
}

exons_sig <- annotate_and_filter(com, com, region_type)
exons_non_sig <- annotate_and_filter(non_sig_sampled, non_sig_sampled, region_type)

# ================================
# Prepare dataframes
# ================================
prepare_df <- function(df, add_group = NULL) {
  df %>%
    mutate(chr = paste0("chr", chr)) %>%
    arrange(chr, start) %>%
    group_by(chr) %>%
    mutate(distance_to_next = lead(start) - end) %>%
    ungroup() %>%
    {if(!is.null(add_group)) mutate(., group = add_group) else .}
}

sig_df <- prepare_df(exons_sig, add_group = "CoM")
non_sig_df <- prepare_df(exons_non_sig, add_group = "No correlation")

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
  xlim(0, 20000) +
  ylim(0, 300)


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
# write.table(close_sig_windows, "genbody_COM_neg_close10k.txt", sep = "\t", row.names = FALSE)
# write.table(close_non_sig_count, "intergenic_no_correlated_distance.txt", sep = "\t", row.names = FALSE)
