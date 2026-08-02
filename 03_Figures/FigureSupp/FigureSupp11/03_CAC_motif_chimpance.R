# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome.Ptroglodytes.NCBI.NHGRImPanTro3v2.0pri)
})

# Define the BSgenome object to avoid repeatedly typing the long name
bsgenome_chimp <- BSgenome.Ptroglodytes.NCBI.NHGRImPanTro3v2.0pri

# --- 1. Load Data and Extract Consensus Regions ---
files <- paste0("~/Escritorio/Paper_scripts/Figure2/Figura2a/data/M", 1:6, "_chimpace_CoM_Windows")

# Unified function to handle both positive (Q3) and negative (Q1) cohorts
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
  
  # Intersect regions across all 6 biological samples
  common_regions <- Reduce(intersect, lapply(processed_list, `[[`, "region_id"))
  
  # Extract the complete rows for the consensus regions (using sample 1 as reference)
  consensus_df <- processed_list[[1]] %>% filter(region_id %in% common_regions)
  return(consensus_df)
}

sig_regions_COM_pos <- extract_consensus_regions(files, "pos")
sig_regions_COM_neg <- extract_consensus_regions(files, "neg")

# Subsample to the lowest common denominator to perfectly prevent class imbalance
min_regions <- min(nrow(sig_regions_COM_pos), nrow(sig_regions_COM_neg))
set.seed(42)
sig_regions_COM_pos <- sig_regions_COM_pos %>% sample_n(min_regions)
sig_regions_COM_neg <- sig_regions_COM_neg %>% sample_n(min_regions)

# --- 2. Chromosome Mapping Dictionary ---
refseq_to_chr <- c(
  "NC_072398.2" = "chr1_hap1_hsa1",   "NC_086015.1" = "chr2_hap1_hsa3",
  "NC_072401.2" = "chr3_hap1_hsa4",   "NC_072402.2" = "chr4_hap1_hsa5",
  "NC_072403.2" = "chr5_hap1_hsa6",   "NC_072404.2" = "chr6_hap1_hsa7",
  "NC_072405.2" = "chr7_hap1_hsa8",   "NC_072406.2" = "chr8_hap1_hsa10",
  "NC_072407.2" = "chr9_hap1_hsa11",  "NC_072408.2" = "chr10_hap1_hsa12",
  "NC_072409.2" = "chr11_hap1_hsa9",  "NC_072410.2" = "chr12_hap1_hsa2a",
  "NC_072411.2" = "chr13_hap1_hsa2b", "NC_072412.2" = "chr14_hap1_hsa13",
  "NC_072413.2" = "chr15_hap1_hsa14", "NC_072414.2" = "chr16_hap1_hsa15",
  "NC_072415.2" = "chr17_hap1_hsa18", "NC_072416.2" = "chr18_hap1_hsa16",
  "NC_072417.2" = "chr19_hap1_hsa17", "NC_072418.2" = "chr20_hap1_hsa19",
  "NC_072419.2" = "chr21_hap1_hsa20", "NC_072420.2" = "chr22_hap1_hsa21",
  "NC_086016.1" = "chr23_hap1_hsa22", "NC_072421.2" = "chrX_hap1_hsaX",
  "NC_072422.2" = "chrY_hap2_hsaY"
)

# Apply mapping
sig_regions_COM_pos <- sig_regions_COM_pos %>% mutate(chr = recode(chr, !!!refseq_to_chr))
sig_regions_COM_neg <- sig_regions_COM_neg %>% mutate(chr = recode(chr, !!!refseq_to_chr))

# --- 3. Motif Discovery Function ---
process_regions <- function(df, group_name) {
  # Boundary safety check
  chr_sizes <- seqlengths(bsgenome_chimp)
  df <- df %>% filter(chr %in% names(chr_sizes), end <= chr_sizes[chr])
  
  gr <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
  
  seqs <- getSeq(bsgenome_chimp, gr)
  
  df_out <- as.data.frame(gr) %>%
    mutate(n_CAC = vcountPattern("CAC", seqs), group = group_name)
  
  return(df_out)
}

# Process target groups
df_com_plus <- process_regions(sig_regions_COM_pos, "mCA(mCG)")
df_com_neg  <- process_regions(sig_regions_COM_neg, "umCA(mCG)")

# --- 4. Generate Random Background ---
# Use the exact chromosomes we defined in the dictionary to build the random background
target_chrs <- unname(refseq_to_chr)
chromsizes <- data.frame(
  chr = target_chrs,
  size = seqlengths(bsgenome_chimp)[target_chrs]
) %>% drop_na()

set.seed(42)
n_random <- nrow(sig_regions_COM_neg)

# Sample chromosomes proportionally to their true genomic size
random_chrs <- sample(chromsizes$chr, n_random, replace = TRUE, 
                      prob = chromsizes$size / sum(chromsizes$size))

# Memory-efficient random start generation ensuring 100bp window boundaries
random_starts <- sapply(random_chrs, function(c) {
  max_start <- chromsizes$size[chromsizes$chr == c] - 100
  sample(0:((max_start - 1) %/% 100), 1) * 100 + 1
})

random_gr <- GRanges(seqnames = random_chrs, ranges = IRanges(start = random_starts, width = 100))
random_df <- as.data.frame(random_gr) %>% rename(chr = seqnames)

df_random <- process_regions(random_df, "Random")

# --- 5. Statistical Analysis ---
df_all_box <- bind_rows(df_com_plus, df_com_neg, df_random) %>%
  mutate(group = factor(group, levels = c("mCA(mCG)", "umCA(mCG)", "Random")))

anova_res <- aov(n_CAC ~ group, data = df_all_box)
cat("\n--- ANOVA Results ---\n")
print(summary(anova_res))

tukey_res <- TukeyHSD(anova_res)
cat("\n--- Tukey HSD Results ---\n")
print(tukey_res)

# --- 6. Categorical Proportion Plotting (0, 1, ≥2) ---
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