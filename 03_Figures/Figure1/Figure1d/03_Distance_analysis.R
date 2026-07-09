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
com <- read_delim("Escritorio/Resultados/CHH/Corr_age/COM_pos_filtrated",
                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

chromsizes <- read.table("/home/epi-lab/Escritorio/Database/hg38.chrom.sizes", header=FALSE) # Obtained from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
colnames(chromsizes) <- c("chr", "size")
chromsizes <- chromsizes[grepl("^chr[0-9XY]+$", chromsizes$chr), ]

gene_features <- readTranscriptFeatures("/home/epi-lab/Escritorio/Database/hg38.bed", remove.unusual = TRUE) 
# Bed file obtained from https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=3856151035_MYca0IwU34JARJDrh0lGVmCABp39. Selecting hg38 genome, group: genes and gene predictions, Track: NCBI Refseq all

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
  scale_fill_manual(values = c("CoM" = "#0793e4", "No correlation" = "#babfba")) +
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
# Function to get close windows (<10kb) without duplication
# ================================
get_close_windows <- function(df, max_dist = 10000) {
  gr <- GRanges(df$chr, IRanges(df$start, df$start + 99))
  
  # Find all windows within max_dist
  hits <- findOverlaps(gr, gr, maxgap = max_dist, type = "any")
  
  # Remove self-overlaps
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  
  # Extract unique indices of close windows
  close_indices <- unique(c(queryHits(hits), subjectHits(hits)))
  
  # Return the filtered dataframe without duplicates
  df[close_indices, ]
}

# ================================
# Apply the function to both groups
# ================================
close_sig_windows <- get_close_windows(sig_df, max_dist = 10002)
close_non_sig_windows <- get_close_windows(non_sig_df, max_dist = 10002)

# ================================
# Fisher's test: proportion of close windows (<10kb)
# ================================

# Add a unique ID to each window
sig_df <- sig_df %>% mutate(id = paste0("sig_", row_number()))
non_sig_df <- non_sig_df %>% mutate(id = paste0("non_", row_number()))
df_all <- bind_rows(sig_df, non_sig_df)

# Detect close windows within each group
close_sig <- get_close_windows(sig_df, max_dist = 10000)
close_non_sig <- get_close_windows(non_sig_df, max_dist = 10000)

# Mark windows that belong to close clusters
df_all <- df_all %>%
  mutate(close = id %in% c(close_sig$id, close_non_sig$id))

# Create a 2x2 contingency table: group vs closeness
contingency <- table(df_all$group, df_all$close)
print(contingency)

# Run Fisher's Exact Test
fisher_result_exact <- fisher.test(contingency[2:1, ], alternative = "greater")
print(fisher_result_exact$p.value, digits=20)

gr_windows <- makeGRangesFromDataFrame(close_sig_windows, 
                                       keep.extra.columns = TRUE, 
                                       seqnames.field = "chr", 
                                       start.field = "start", 
                                       end.field = "end")

gr_sig <- makeGRangesFromDataFrame(sig_df, 
                                   seqnames.field = "chr", 
                                   start.field = "start", 
                                   end.field = "end")


gr_far_windows <- gr_sig[!gr_sig %over% gr_windows, ]


close_far_windows <- as.data.frame(gr_far_windows)
# Optionally save results
write.table(close_sig_windows, "promoters_COM_neg_close10k.txt", sep = "\t", row.names = FALSE)
write.table(close_far_windows, "promoters_COM_neg_far10k.txt", sep = "\t", row.names = FALSE)

#############################################################################################
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(readr)

# ---------------------------
# Load CoM windows
# ---------------------------
com_pos <- read_delim(
  "Escritorio/Resultados/CHH/Corr_age/COM_pos_filtrated",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
) %>% mutate(chr = paste0("chr", chr))

com_neg <- read_delim(
  "Escritorio/Resultados/CHH/Corr_age/COM_neg_filtrated",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
) %>% mutate(chr = paste0("chr", chr))

# ---------------------------
# Convert to GRanges
# ---------------------------
gr_pos <- makeGRangesFromDataFrame(com_pos, keep.extra.columns = TRUE)
gr_neg <- makeGRangesFromDataFrame(com_neg, keep.extra.columns = TRUE)

# ---------------------------
# Distance from negative windows to nearest positive window
# ---------------------------
nearest_pos <- nearest(gr_neg, gr_pos)  
dist_neg_to_pos <- abs(start(gr_pos[nearest_pos]) - end(gr_neg))

# ---------------------------
# Generate random windows for background
# ---------------------------
set.seed(42)
chrom_sizes <- read.table("/home/epi-lab/Escritorio/Database/hg38.chrom.sizes", header=FALSE)
colnames(chrom_sizes) <- c("chr", "size")
chrom_sizes <- chrom_sizes[grepl("^chr[0-9XY]+$", chrom_sizes$chr), ]

# Sample random chromosomes proportional to size
random_chrs <- sample(
  chrom_sizes$chr,
  length(gr_neg),
  replace = TRUE,
  prob = chrom_sizes$size / sum(chrom_sizes$size)
)

# Generate random starts
chrom_start_list <- setNames(
  lapply(chrom_sizes$size, function(sz) seq(1, sz - 100, by = 100)),
  chrom_sizes$chr
)

random_starts <- mapply(
  function(chr, n) sample(chrom_start_list[[chr]], n, replace = TRUE),
  chr = random_chrs, n = 1
)

random_gr <- GRanges(seqnames = random_chrs, ranges = IRanges(start=random_starts, width=100))

# Distance from negative windows to random windows
nearest_random <- nearest(gr_neg, random_gr)
dist_neg_to_random <- abs(start(random_gr[nearest_random]) - end(gr_neg))

# ---------------------------
# Combine distances
# ---------------------------
df_dist <- data.frame(
  distance = c(dist_neg_to_pos, dist_neg_to_random),
  group = rep(c("CoM", "Not Correlated"),
              times = c(length(dist_neg_to_pos), length(dist_neg_to_random)))
)

# ---------------------------
# Plot histogram
# ---------------------------
ggplot(df_dist, aes(x = distance, fill = group)) +
  geom_histogram(binwidth = 100, position = "identity", alpha = 1) +
  scale_fill_manual(values = c("CoM" = "#f28229", "Not Correlated" = "#babfba"))  +
  labs(
    x = "Distance (bp)",
    y = "Frecuency",
    fill = NULL
  )+
  theme_linedraw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 1.5, colour = "black"),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.ticks.y = element_line(size = 1),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_line(size = 1),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) +
  xlim(0, 20000)+ylim(0,300)

# ====================================================
# Fisher exact test: proporción de ventanas cercanas
# ====================================================

# Definir umbral de distancia como "cercano"
threshold <- 10000  # puedes ajustar si quieres

df_dist$close <- df_dist$distance <= threshold

# Crear la tabla 2x2
contingency <- table(df_dist$group, df_dist$close)
print("Contingency table:")
print(contingency)

# Aseguramos orden consistente
df_dist$group <- factor(df_dist$group, levels = c("CoM", "Not Correlated"))
contingency <- table(df_dist$group, df_dist$close)

# Fisher exact test
fisher_res <- fisher.test(contingency, alternative = "greater")

cat("\n=== Fisher Exact Test Results ===\n")
cat("Odds Ratio:", fisher_res$estimate, "\n")
cat("IC95%:", fisher_res$conf.int, "\n")
cat("p-value:", fisher_res$p.value, "\n\n")
