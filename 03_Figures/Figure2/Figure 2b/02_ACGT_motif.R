
# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(tidyverse)
  
})

# Define input file paths
PATH_COM_POS <- "/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/genbody_COM_pos_close10k.txt"
PATH_COM_NEG <- "/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/genbody_COM_neg_close10k.txt"
path_universe_windows   <- "~/Escritorio/Resultados/CHH/Corr_age/all_regions.txt"
# --- 1. Data Ingestion & Preprocessing ---
# Load significant regions correlated with age (CoM+ and CoM-)
sig_regions_COM_pos <- read_delim(PATH_COM_POS, delim = "\t", show_col_types = FALSE)
sig_regions_COM_neg <- read_delim(PATH_COM_NEG, delim = "\t", show_col_types = FALSE)
universe_windows_df <- read_delim(path_universe_windows, delim = "\t", trim_ws = TRUE)
universe_windows_df$chr = paste0("chr", universe_windows_df$chr)
# Subsample CoM+ to match the exact size of CoM- to prevent class imbalance biases

com_pos_q3 <- sig_regions_COM_pos %>% filter(diff > quantile(diff, 0.75, na.rm = TRUE))
com_neg_q3 <- sig_regions_COM_neg %>% filter(diff < quantile(diff, 0.75, na.rm = TRUE))
set.seed(42)
com_neg_q3 <- com_neg_q3 %>% sample_n(size = nrow(com_pos_q3), replace = FALSE)


# --- 2. Motif Discovery (ACTG) ---
# Function to count absolute occurrences of 'ACTG' motifs in canonical chromosomes
process_regions <- function(df, group_name) {
  gr <- GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start, end = df$end)
  )
  
  # Restrict analysis to canonical chromosomes (1-22, X, Y, M)
  canonical_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
  gr <- gr[seqnames(gr) %in% canonical_chrs]
  
  # Extract sequences and count ACTG motifs
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  n_ACTG <- vcountPattern("ACTG", seqs)
  
  # Format output
  df_out <- as.data.frame(gr) %>%
    mutate(
      n_ACTG = n_ACTG,
      group = group_name
    )
  return(df_out)
}

# --- 4. Process Target Groups and Generate Null Background ---
df_com_plus <- process_regions(com_pos_q3, "CoM +")
df_com_neg  <- process_regions(com_neg_q3, "CoM -")

# Generate size-matched background (Not Correlated) from the overarching universe dataset
set.seed(42)
df_random <- universe_windows_df %>% sample_n(size = nrow(com_neg_q3), replace = FALSE)
df_random <- process_regions(df_random, "Not Correlated")

# --- 5. Statistical Analysis ---
# Combine all cohorts into a single dataframe for statistical modeling
df_all_box <- bind_rows(df_com_plus, df_com_neg, df_random) %>%
  mutate(group = factor(group, levels = c("CoM +", "CoM -", "Not Correlated")))

# Perform One-Way ANOVA to test for variance in ACTG motifs across groups
anova_res <- aov(n_ACTG ~ group, data = df_all_box)
cat("\n--- ANOVA Results ---\n")
print(summary(anova_res))

# Perform Tukey's HSD for pairwise group comparisons
tukey_res <- TukeyHSD(anova_res)
cat("\n--- Tukey HSD Results ---\n")
print(tukey_res)

# --- 6. Plotting Preparation ---
# Categorize motif counts into discrete bins for visualization
df_filtered <- df_all_box %>%
  mutate(
    cat_plot = case_when(
      n_ACTG == 0 ~ "0 motifs",
      n_ACTG == 1 ~ "1 motif",
      n_ACTG >= 2 ~ "≥2 motifs"
    ),
    cat_plot = factor(cat_plot, levels = c("0 motifs", "1 motif", "≥2 motifs"))
  )

# Calculate relative proportions of motif categories within each group
df_prop_plot <- df_filtered %>%
  group_by(cat_plot, group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cat_plot) %>% # Standardizing proportions relative to the motif category
  mutate(prop_within_cat = count / sum(count)) %>%
  ungroup()

# Calculate specific risk/enrichment ratio (≥2 motifs vs 0 motifs) for CoM+
df_coMplus <- df_prop_plot %>%
  filter(group == "CoM -", cat_plot %in% c("0 motifs", "≥2 motifs"))

prop_2 <- df_coMplus$prop_within_cat[df_coMplus$cat_plot == "≥2 motifs"]
prop_0 <- df_coMplus$prop_within_cat[df_coMplus$cat_plot == "0 motifs"]

ratio_2_vs_0 <- prop_2 / prop_0
cat("\nRatio (≥2 motifs / 0 motifs) for CoM -:", ratio_2_vs_0, "\n")

# --- 7. Visualization ---
# Ensure correct factor leveling for aesthetics
df_prop_plot$group <- factor(df_prop_plot$group, levels = c("CoM -", "CoM +", "Not Correlated"))

publication_plot <- ggplot(df_prop_plot, aes(x = cat_plot, y = prop_within_cat, fill = group)) +
  geom_col(alpha = 0.7)+
  scale_fill_manual(values = c("CoM +" = "#0072B2", "CoM -" = "#C22506", "Not Correlated" = "#7c807c")) +
  labs(
    x = "Nº of ACTG motif",
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

# Display the plot
print(publication_plot)
