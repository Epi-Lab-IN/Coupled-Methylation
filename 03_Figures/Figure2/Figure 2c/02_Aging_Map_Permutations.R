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
  attributes = c('ensembl_gene_id','chromosome_name','start_position','end_position','gene_biotype'),
  mart = ensembl
)
standard_chrs <- c(as.character(1:22), "X", "Y", "MT")
all_genes <- all_genes[all_genes$chromosome_name %in% standard_chrs, ]
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
# 3. LOAD SIGNIFICANT REGIONS & SUBSAMPLING
# ================================================================
# Load data
com_pos <- read_delim("/home/epi-lab/Escritorio/Resultados/CHH/genes/genes_COM_pos_aggregated.txt", delim = "\t", show_col_types = FALSE)
com_neg <- read_delim("/home/epi-lab/Escritorio/Resultados/CHH/genes/genes_COM_neg_aggregated.txt", delim = "\t", show_col_types = FALSE)

# Annotate and filter genes with >= 2 regions
genes_pos <- com_pos %>% dplyr::filter(n_regions >= 2)
genes_neg <- com_neg %>% dplyr::filter(n_regions >= 2)

min_genes_count <- nrow(genes_pos) 

# ================================================================
# 4. MAP TO AGINGMAP
# ================================================================
aging_map <- read_delim("~/Escritorio/Database/AgingMap_Humans.txt", delim="\t", trim_ws=TRUE, show_col_types = FALSE) #Select database for compare (https://sysomics.com/AgingMap/,https://epifactors.autosome.org, Kegg_AD )

# Merge target genes with HGNC symbols
annot_sig <- merge(genes_pos, annot_all, by.x="ensembl_gene_id", by.y="input")

# Calculate observed overlap with AgingMap database
sig_genes_in_aging <- merge(annot_sig, aging_map, by.x="name", by.y="Symbol")
observed_overlap <- nrow(sig_genes_in_aging)

cat("Observed overlap with AgingMap:", observed_overlap, "\n")
# ================================================================
# 5. PERMUTATIONS (direct random gene sampling)
# ================================================================
n_perm <- 1000
set.seed(455)
overlap_counts <- numeric(n_perm)

pb <- txtProgressBar(min=0, max=n_perm, style=3)

for(i in 1:n_perm){
  # 1. Randomly sample the same number of genes from the full gene set
  rand_ids <- sample(all_gene_ids, min_genes_count, replace = FALSE)
  
  # 2. Get HGNC symbols for the sampled genes
  annot_rand <- annot_all %>% dplyr::filter(input %in% rand_ids)
  
  # 3. Count overlaps with AgingMap
  overlap_counts[i] <- nrow(merge(annot_rand, aging_map, by.x="name", by.y="Symbol"))
  
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

