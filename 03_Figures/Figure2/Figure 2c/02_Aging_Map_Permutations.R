# ================================================================
# LIBRARIES & SETUP
# ================================================================
suppressPackageStartupMessages({
  library(readr)
  library(biomaRt)
  library(GenomicRanges)
  library(dplyr)
  library(gprofiler2)
})

# ================================================================
# 1. GET GENE ANNOTATION (ONE TIME)
# ================================================================
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
all_genes <- getBM(
  attributes = c('ensembl_gene_id','chromosome_name','start_position','end_position'),
  filters = 'biotype', 
  values = "protein_coding", 
  mart = ensembl
)

all_genes$chromosome_name <- paste0("chr", all_genes$chromosome_name)
genes_gr <- makeGRangesFromDataFrame(
  all_genes,
  seqnames.field = "chromosome_name",
  start.field = "start_position",
  end.field = "end_position",
  keep.extra.columns = TRUE
)

# Extract all genes present in the annotation for stable ID conversion
all_gene_ids <- unique(genes_gr$ensembl_gene_id)
annot_all <- gconvert(all_gene_ids, organism="hsapiens")

# ================================================================
# 2. OPTIMIZED ANNOTATION FUNCTION
# ================================================================

annotate_regions <- function(regions_df, genes_gr, mode = c("min", "max", "random")) {
  mode <- match.arg(mode)
  if(nrow(regions_df) == 0) return(data.frame())
  
  regions_gr <- makeGRangesFromDataFrame(regions_df, keep.extra.columns = TRUE)
  hits <- findOverlaps(regions_gr, genes_gr)
  if(length(hits) == 0) return(data.frame())
  
  
  df <- data.frame(
    ensembl_gene_id = genes_gr$ensembl_gene_id[subjectHits(hits)]
  )
  
  # Behavior based on the selected mode
  if (mode == "random") {
    res <- df %>%
      group_by(ensembl_gene_id) %>%
      summarise(n_regions = n(), .groups = "drop")
    return(as.data.frame(res))
  } 
  
  df$diff <- regions_gr$diff[queryHits(hits)]
  
  res <- df %>%
    group_by(ensembl_gene_id) %>%
    summarise(
      target_diff = if(mode == "min") min(diff) else max(diff),
      n_regions = n(),
      .groups = "drop"
    ) %>%
    arrange(if(mode == "min") target_diff else dplyr::desc(target_diff))
  
  return(as.data.frame(res))
}

# ================================================================
# 3. LOAD SIGNIFICANT REGIONS & SUBSAMPLING
# ================================================================
# Load data
com_pos <- read_delim("/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/genbody_COM_pos_close10k.txt", delim = "\t", show_col_types = FALSE)
com_neg <- read_delim("/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/genbody_COM_neg_close10k.txt", delim = "\t", show_col_types = FALSE)

# Annotate and filter genes with >= 2 regions
genes_pos <- annotate_regions(com_pos, genes_gr, mode = "max") %>% dplyr::filter(n_regions >= 2)
genes_neg <- annotate_regions(com_neg, genes_gr, mode = "min") %>% dplyr::filter(n_regions >= 2)

min_genes_count <- nrow(genes_neg) 

# ================================================================
# 4. MAP TO AGINGMAP
# ================================================================
aging_map <- read_delim("~/Escritorio/Database/epigenetic_genes.txt", delim="\t", trim_ws=TRUE, show_col_types = FALSE) #Select database for compare (https://sysomics.com/AgingMap/,https://epifactors.autosome.org, Kegg_AD )

# Merge target genes with HGNC symbols
annot_sig <- merge(genes_neg, annot_all, by.x="ensembl_gene_id", by.y="input")

# Calculate observed overlap with AgingMap database
sig_genes_in_aging <- merge(annot_sig, aging_map, by.x="name", by.y="Symbol")
observed_overlap <- nrow(sig_genes_in_aging)

cat("Observed overlap with AgingMap:", observed_overlap, "\n")
# ================================================================
# 5. PERMUTATIONS
# ================================================================
n_perm <- 1000
set.seed(42)
overlap_counts <- numeric(n_perm)



chrom_sizes <- read.table("/home/epi-lab/Escritorio/Database/hg38.chrom.sizes", header=FALSE)
colnames(chrom_sizes) <- c("chr", "size")
chrom_sizes <- chrom_sizes[grepl("^chr[0-9XY]+$", chrom_sizes$chr), ]
pb <- txtProgressBar(min=0, max=n_perm, style=3)

for(i in 1:n_perm){
  # 1. Select random chromosomes weighted by their size
  rand_chrs <- sample(chrom_sizes$chr, min_genes_count, replace=TRUE, 
                      prob=chrom_sizes$size / sum(chrom_sizes$size))
  
  # 2. Select random start positions
  rand_starts <- sapply(rand_chrs, function(chr) {
    max_start <- chrom_sizes$size[chrom_sizes$chr == chr] - 100
    sample(1:max_start, 1)
  })
  
  # 3. Create random GRanges object
  rand_df <- data.frame(chr = rand_chrs, start = rand_starts, end = rand_starts + 100)
  
  # 4. Annotate random regions and filter by >= 2 regions
  rand_genes <- annotate_regions(rand_df, genes_gr, mode = "random") %>% 
    filter(n_regions >= 2)
  
  # 5. Cap to min_genes_count if exceeded; otherwise use all available genes
  if(nrow(rand_genes) > min_genes_count){
    rand_genes <- rand_genes %>% sample_n(min_genes_count)
  }
  
  # 6. Count overlaps with AgingMap
  if(nrow(rand_genes) > 0){
    annot_rand <- merge(rand_genes, annot_all, by.x="ensembl_gene_id", by.y="input")
    overlap_counts[i] <- nrow(merge(annot_rand, aging_map, by.x="name", by.y="Symbol"))
  } else {
    overlap_counts[i] <- 0
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

# ================================================================
# 6. STATISTICS
# ================================================================
empirical_p <- sum(overlap_counts >= observed_overlap) / n_perm
enrichment <- observed_overlap / mean(overlap_counts)

cat("\n--- RESULTS ---\n")
cat("Empirical p-value:", empirical_p, "\n")
cat("Enrichment factor:", round(enrichment, 2), "\n")