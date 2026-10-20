# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome.Drerio.NCBI.GRCz12tu)
})

# Define the BSgenome object
bsgenome_zebrafish <- BSgenome.Drerio.NCBI.GRCz12tu

# Define input file paths 
PATH_WINDOWS <- "~/Escritorio/Paper_scripts/Figure2/Figura2c/Data/Zebrafish/Zebrafish_CoM_Windows"
PATH_MAPPING <- "~/Escritorio/Animals/Genomes/Assembly_zebra.txt"

# --- 1. Load Data and Format Chromosomes ---
sample_data <- read_delim(PATH_WINDOWS, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
map_data <- read_delim(PATH_MAPPING, delim = "\t", show_col_types = FALSE)

# Create mapping dictionary (RefSeq <-> Scaffold)
dicc <- map_data %>%
  dplyr::select(refseq = GenBankseq, scaffold = Chromosomename)

# Apply mapping, format with "chr" prefix, and drop unmapped regions safely
sample_data <- sample_data %>%
  left_join(dicc, by = c("chr" = "refseq")) %>%
  mutate(chr = paste0("chr", scaffold)) %>%
  dplyr::select(-scaffold) %>%
  filter(!is.na(chr) & chr != "chrNA")

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
  gr <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
  
  
  chr_lengths <- seqlengths(bsgenome_zebrafish)
  valid_chrs <- names(chr_lengths)
  
  # Filter for valid chromosomes and boundaries
  gr <- gr[seqnames(gr) %in% valid_chrs]
  gr <- gr[end(gr) <= chr_lengths[as.character(seqnames(gr))] & start(gr) >= 1]
  gr <- trim(gr, use.names = TRUE)
  
  # Extract sequences and count absolute occurrences of 'CAC'
  seqs <- getSeq(bsgenome_zebrafish, gr)
  
  df_out <- as.data.frame(gr) %>%
    mutate(n_CAC = vcountPattern("CAC", seqs), group = group_name)
  
  return(df_out)
}

# Process target groups
df_com_plus <- process_regions(sig_regions_COM_pos, "mCA(mCG)")
df_com_neg  <- process_regions(sig_regions_COM_neg, "umCA(mCG)")

# --- 4. Generate Random Background ---
# Filter strictly for main canonical chromosomes (chr1 to chr25)
bs_chrnames <- seqnames(bsgenome_zebrafish)
bs_chrnames_main <- bs_chrnames[bs_chrnames %in% paste0("chr", 1:25)]

chromsizes <- data.frame(
  chr = bs_chrnames_main,
  size = seqlengths(bsgenome_zebrafish)[bs_chrnames_main]
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
    axis.line = element_line(size = 0.25),  # Grosor de los ejes
    axis.ticks = element_line(size = 0.25) # Grosor de las marcas
  )