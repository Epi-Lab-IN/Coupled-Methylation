
# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome.Xtropicalis.NCBI.UCBXtro10.0)
  # library(biomaRt) # Unused in this specific script, removed to save memory
})

# Define the BSgenome object
bsgenome_frog <- BSgenome.Xtropicalis.NCBI.UCBXtro10.0

# Define input file paths 
PATH_WINDOWS <- "/home/epi-lab/Escritorio/Paper_scripts/Figure2/Figura2c/Data/Frog/Frog_CoM_Windows"

# --- 1. Load Data and Format Chromosomes ---
sample_data <- read_delim(PATH_WINDOWS, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# Prepend "Chr" to match the BSgenome formatting for X. tropicalis
sample_data <- sample_data %>% mutate(chr = paste0("Chr", chr))

# --- 2. Filter Significant Regions (CoM+ and CoM-) ---
# Calculate quantiles for CA methylation values
q3 <- quantile(sample_data$CA_bvalues, 0.75, na.rm = TRUE)
q1 <- quantile(sample_data$CA_bvalues, 0.25, na.rm = TRUE)

# Isolate target regions based on methylation thresholds and CG > 66
sig_regions_COM_pos <- sample_data %>% filter(CG_bvalues > 66 & CA_bvalues > q3)
sig_regions_COM_neg <- sample_data %>% filter(CG_bvalues > 66 & CA_bvalues <= q1)

# Safely subsample to the lowest common denominator to perfectly prevent class imbalance
min_regions <- min(nrow(sig_regions_COM_pos), nrow(sig_regions_COM_neg))

set.seed(42)
sig_regions_COM_pos <- sig_regions_COM_pos %>% sample_n(min_regions)
sig_regions_COM_neg <- sig_regions_COM_neg %>% sample_n(min_regions)

# --- 3. Motif Discovery Function ---
process_regions <- function(df, group_name) {
  # Boundary safety check: ensure regions don't spill off the end of a chromosome
  chr_sizes <- seqlengths(bsgenome_frog)
  df <- df %>% filter(chr %in% names(chr_sizes), end <= chr_sizes[chr])
  
  gr <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
  
  # Extract sequences and count absolute occurrences of 'CAC'
  seqs <- getSeq(bsgenome_frog, gr)
  
  df_out <- as.data.frame(gr) %>%
    mutate(n_CAC = vcountPattern("CAC", seqs), group = group_name)
  
  return(df_out)
}

# Process target groups
df_com_plus <- process_regions(sig_regions_COM_pos, "mCA(mCG)")
df_com_neg  <- process_regions(sig_regions_COM_neg, "umCA(mCG)")

# --- 4. Generate Random Background ---
# Filter for main canonical chromosomes 
bs_chrnames <- seqnames(bsgenome_frog)
bs_chrnames_main <- bs_chrnames[grepl("^Chr[0-9]+$", bs_chrnames)]

chromsizes <- data.frame(
  chr = bs_chrnames_main,
  size = seqlengths(bsgenome_frog)[bs_chrnames_main]
)

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

# --- 5. Statistical Analysis ---
df_all_box <- bind_rows(df_com_plus, df_com_neg, df_random) %>%
  mutate(group = factor(group, levels = c("mCA(mCG)", "umCA(mCG)", "Random")))

anova_res <- aov(n_CAC ~ group, data = df_all_box)
cat("\n--- ANOVA Results ---\n")
print(summary(anova_res))

tukey_res <- TukeyHSD(anova_res)
cat("\n--- Tukey HSD Results ---\n")
print(tukey_res)

# --- 6. Categorical Motif Enrichment Plotting (0, 1, ≥2) ---
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
    axis.line = element_line(size = 0.25), 
    axis.ticks = element_line(size = 0.25) 
  )
