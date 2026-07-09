# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# Define the BSgenome object
bsgenome_mouse <- BSgenome.Mmusculus.UCSC.mm10

# Define input file paths 
PATH_WINDOWS_DIR <- "Escritorio/Animals/CHH/Mouse/"
PATH_CHROM_SIZES <- "/home/epi-lab/Escritorio/Database/mm10.chrom.sizes"

files <- paste0(PATH_WINDOWS_DIR, "cowindows_", c("6w", "10w", "JC", "VC"), ".txt")


# --- 1. Load Data and Extract Consensus Regions ---
extract_consensus_regions <- function(file_list, type = c("pos", "neg")) {
  type <- match.arg(type)
  
  processed_list <- lapply(file_list, function(f) {
    df <- read_delim(f, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
    
    if (type == "pos") {
      threshold <- quantile(df$CA_bvalues, 0.75, na.rm = TRUE)
      df_filt <- df %>% filter(CA_bvalues >= threshold & CG_bvalues > 66)
    } else {
      threshold <- quantile(df$CA_bvalues, 0.25, na.rm = TRUE)
      df_filt <- df %>% filter(CA_bvalues <= threshold & CG_bvalues > 66)
    }
    
    # Create unique identifier for regions
    df_filt %>% mutate(region_id = paste(chr, start, end, sep = "_"))
  })
  
  # Intersect regions across all samples (6w, 10w, 24w, 72w)
  common_regions <- Reduce(intersect, lapply(processed_list, `[[`, "region_id"))
  
  # Extract the complete rows for the consensus regions
  consensus_df <- processed_list[[1]] %>% filter(region_id %in% common_regions)
  return(consensus_df)
}

sig_regions_COM_pos <- extract_consensus_regions(files, "pos")
sig_regions_COM_neg <- extract_consensus_regions(files, "neg")

# Safely subsample to the lowest common denominator to perfectly prevent class imbalance
min_regions <- min(nrow(sig_regions_COM_pos), nrow(sig_regions_COM_neg))

set.seed(42)
sig_regions_COM_pos <- sig_regions_COM_pos %>% sample_n(min_regions)
sig_regions_COM_neg <- sig_regions_COM_neg %>% sample_n(min_regions)

# --- 2. Motif Discovery Function ---
process_regions <- function(df, group_name) {
  gr <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
  
  # Filter strictly for main mouse autosomes (chr1 to chr19)
  canonical_chrs <- paste0("chr", 1:19)
  gr <- gr[seqnames(gr) %in% canonical_chrs]
  
  # Extract sequences and count absolute occurrences of 'CAC'
  seqs <- getSeq(bsgenome_mouse, gr)
  
  df_out <- as.data.frame(gr) %>%
    mutate(n_CAC = vcountPattern("CAC", seqs), group = group_name)
  
  return(df_out)
}

# Process target groups
df_com_plus <- process_regions(sig_regions_COM_pos, "mCA(mCG)")
df_com_neg  <- process_regions(sig_regions_COM_neg, "umCA(mCG)")

# --- 3. Generate Random Background ---
chromsizes <- read.table(PATH_CHROM_SIZES, header = FALSE, col.names = c("chr", "size"))

# Ensure the random background ONLY uses the canonical autosomes (chr1-19) to match real data
canonical_chrs <- paste0("chr", 1:19)
chromsizes <- chromsizes %>% filter(chr %in% canonical_chrs)

set.seed(42)
n_random <- min_regions

# Sample chromosomes proportionally to their true genomic size
random_chrs <- sample(chromsizes$chr, n_random, replace = TRUE, 
                      prob = chromsizes$size / sum(chromsizes$size))

# Memory-efficient random start generation ensuring 100bp window boundaries
random_starts <- sapply(random_chrs, function(c) {
  max_start <- chromsizes$size[chromsizes$chr == c] - 100
  sample(0:((max_start - 1) %/% 100), 1) * 100 + 1
})

# Construct GRanges and process motifs
random_gr <- GRanges(seqnames = random_chrs, ranges = IRanges(start = random_starts, width = 100))
random_df <- as.data.frame(random_gr) %>% rename(chr = seqnames)

df_random <- process_regions(random_df, "Random")

# --- 4. Statistical Analysis ---
df_all_box <- bind_rows(df_com_plus, df_com_neg, df_random) %>%
  mutate(group = factor(group, levels = c("mCA(mCG)", "umCA(mCG)", "Random")))

anova_res <- aov(n_CAC ~ group, data = df_all_box)
cat("\n--- ANOVA Results ---\n")
print(summary(anova_res))

tukey_res <- TukeyHSD(anova_res)
cat("\n--- Tukey HSD Results ---\n")
print(tukey_res)

# --- 5. Categorical Motif Enrichment Plotting (0, 1, ≥2) ---
df_filtered <- df_all_box %>%
  mutate(
    cat_plot = case_when(
      n_CAC == 0 ~ "0 motifs",
      n_CAC == 1 ~ "1 motif",
      n_CAC >= 2 ~ "≥2 motifs"
    ),
    cat_plot = factor(cat_plot, levels = c("0 motifs", "1 motif", "≥2 motifs"))
  )

# Calculate proportions by group and category
df_prop_plot <- df_filtered %>%
  group_by(cat_plot, group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cat_plot) %>%
  mutate(prop_within_cat = count / sum(count)) %>%
  ungroup()

# Calculate specific ratio (≥2 / 0) for mCA(mCG)
df_coMplus <- df_prop_plot %>% filter(group == "mCA(mCG)")

prop_2 <- df_coMplus$prop_within_cat[df_coMplus$cat_plot == "≥2 motifs"]
prop_0 <- df_coMplus$prop_within_cat[df_coMplus$cat_plot == "0 motifs"]

ratio_2_vs_0 <- prop_2 / prop_0
cat("\nRatio (≥2 motifs / 0 motifs) for mCA(mCG):", ratio_2_vs_0, "\n")

# Generate the plot
ggplot(df_prop_plot, aes(x = cat_plot, y = prop_within_cat, fill = group)) +
  geom_col(alpha = 0.7)+
  scale_fill_manual(values = c("mCA(mCG)" = "#0072B2", "umCA(mCG)" = "#C22506", "Random" = "#7c807c")) +
  labs(
    x = "Nº of CAC motif",
    y = "Proportion"
  ) +
  theme_linedraw(base_size = 14)+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.ticks.y = element_line(size = 1),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_line(size = 1),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    axis.line = element_line(size = 0.25),  # Grosor de los ejes
    axis.ticks = element_line(size = 0.25) # Grosor de las marcas
  )