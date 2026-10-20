
# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome.Ggallus.NCBI.bGalGal1.mat.broiler.GRCg7b)
  # library(biomaRt) # Unused in this specific script, removed to save memory
})

# Define input file paths (Update when publishing)
PATH_WINDOWS <- "/home/epi-lab/Escritorio/Paper_scripts/Figure2/Figura2c/Data/Chicken/Chicken_CoM_Windows"

# --- 1. Load and Map Genomic Data ---
sample_data <- read_delim(PATH_WINDOWS, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# Define main contigs
contigs_main <- c("CM028482.1","CM028483.1","CM028484.1","CM028485.1",
                  "CM028486.1","CM028487.1","CM028488.1","CM028489.1",
                  "CM028490.1","CM028491.1","CM028492.1","CM028493.1",
                  "CM028494.1","CM028495.1","CM028496.1","CM028497.1",
                  "CM028498.1","CM028499.1","CM028500.1","CM028501.1",
                  "CM028502.1","CM028503.1","CM028504.1","CM028505.1",
                  "CM028506.1","CM028507.1","CM028508.1","CM028509.1",
                  "CM028510.1","CM028511.1","CM028512.1","CM028513.1",
                  "CM028514.1","CM028515.1","CM028516.1","CM028517.1",
                  "CM028518.1","CM028519.1","CM028520.1","CM028521.1",
                  "CM028522.1","CM028585.1") 

# Map contigs to BSgenome scaffold names
bs_names <- seqnames(BSgenome.Ggallus.NCBI.bGalGal1.mat.broiler.GRCg7b)
bs_main <- bs_names[1:length(contigs_main)]
map_scaffold <- setNames(bs_main, contigs_main)

# Apply mapping and immediately drop regions with NA chromosomes
sample_data <- sample_data %>%
  mutate(chr = map_scaffold[chr]) %>%
  filter(!is.na(chr))

# --- 2. Filter Significant Regions (CoM+ and CoM-) ---
# Calculate quantiles for CA methylation values
q3 <- quantile(sample_data$CA_bvalues, 0.75, na.rm = TRUE)
q1 <- quantile(sample_data$CA_bvalues, 0.25, na.rm = TRUE)

sig_regions_COM_pos <- sample_data %>% filter(CG_bvalues > 66 & CA_bvalues > q3)
sig_regions_COM_neg <- sample_data %>% filter(CG_bvalues > 66 & CA_bvalues <= q1)

# Subsample CoM+ to match the exact size of CoM- (prevents class imbalance)
set.seed(42)
sig_regions_COM_pos <- sig_regions_COM_pos %>% sample_n(nrow(sig_regions_COM_neg), replace = FALSE)

# --- 3. Motif Discovery Function ---
process_regions <- function(df, group_name) {
  chr_sizes <- seqlengths(BSgenome.Ggallus.NCBI.bGalGal1.mat.broiler.GRCg7b)
  df <- df %>% filter(chr %in% names(chr_sizes), end <= chr_sizes[chr])
  
  gr <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
  
  seqs <- getSeq(BSgenome.Ggallus.NCBI.bGalGal1.mat.broiler.GRCg7b, gr)
  n_CAC <- vcountPattern("CAC", seqs)
  
  df_out <- as.data.frame(gr) %>%
    mutate(n_CAC = n_CAC, group = group_name)
  return(df_out)
}

df_com_plus <- process_regions(sig_regions_COM_pos, "mCA(mCG)")
df_com_neg  <- process_regions(sig_regions_COM_neg, "umCA(mCG)")

# --- 4. Generate Random Background (Memory Optimized) ---
valid_chr <- seqnames(BSgenome.Ggallus.NCBI.bGalGal1.mat.broiler.GRCg7b)
chromsizes <- seqlengths(BSgenome.Ggallus.NCBI.bGalGal1.mat.broiler.GRCg7b)[valid_chr]

chromsizes_df <- data.frame(chr = names(chromsizes), size = as.numeric(chromsizes))

set.seed(42)
n_random <- nrow(sig_regions_COM_neg)

random_chrs <- sample(chromsizes_df$chr, n_random, replace = TRUE, 
                      prob = chromsizes_df$size / sum(chromsizes_df$size))


random_starts <- sapply(random_chrs, function(c) {
  max_start <- chromsizes_df$size[chromsizes_df$chr == c] - 100
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
    # Standardized level ordering from 0 to >=2 for a natural x-axis progression
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

