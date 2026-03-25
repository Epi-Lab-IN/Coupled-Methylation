
# --- 0. Setup and Dependencies ---
suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(tidyverse)
  
})

# Define input file paths
PATH_COM_POS <- "/home/epi-lab/Escritorio/Resultados/Corr_with_age/COM_pos_close10k"
PATH_COM_NEG <- "/home/epi-lab/Escritorio/Resultados/Corr_with_age/COM_neg_close10k"
PATH_CA_METH <- "/home/epi-lab/Escritorio/Resultados/ventanas_ca_comun_cov5_bvalues"

# --- 1. Data Ingestion & Preprocessing ---
# Load significant regions correlated with age (CoM+ and CoM-)
sig_regions_COM_pos <- read_delim(PATH_COM_POS, delim = "\t", show_col_types = FALSE)
sig_regions_COM_neg <- read_delim(PATH_COM_NEG, delim = "\t", show_col_types = FALSE)

# Isolate the top quartiles based on correlation differences
com_pos_q3 <- sig_regions_COM_pos %>% filter(diff > quantile(diff, 0.75, na.rm = TRUE))
com_neg_q3 <- sig_regions_COM_neg %>% filter(diff < quantile(diff, 0.75, na.rm = TRUE))

# Subsample CoM+ to match the exact size of CoM- to prevent class imbalance biases
set.seed(42)
com_pos_q3_subsampled <- com_pos_q3 %>% sample_n(size = nrow(com_neg_q3), replace = FALSE)

# Load non-CpG (CA) methylation data
ca_data <- read_delim(PATH_CA_METH, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# Format chromosome names and calculate mean non-CpG methylation across target columns (assumed indices 4:10)
ca_data <- ca_data %>%
  mutate(
    chr = paste0("chr", chr),
    mean_ca = rowMeans(dplyr::select(., 4:10), na.rm = TRUE)
  )

# --- 2. Genomic Ranges Conversion & Overlap Analysis ---
# Convert CA methylation data to GRanges object
gr_ca <- GRanges(
  seqnames = ca_data$chr,
  ranges   = IRanges(start = ca_data$start, end = ca_data$end),
  mean_ca  = ca_data$mean_ca
)

# Function to map mean CA methylation to significant regions via genomic overlaps
process_hits <- function(df_regions, gr_target) {
  gr_regions <- GRanges(
    seqnames = df_regions$chr,
    ranges   = IRanges(start = df_regions$start, end = df_regions$end)
  )
  
  # Identify overlaps between query regions and target CA regions
  hits <- findOverlaps(gr_regions, gr_target)
  
  # Assign the corresponding mean CA methylation; non-overlapping regions get NA
  df_regions$mean_ca <- NA
  df_regions$mean_ca[queryHits(hits)] <- gr_target$mean_ca[subjectHits(hits)]
  
  return(df_regions)
}

sig_regions_COM_pos <- process_hits(com_pos_q3_subsampled, gr_ca)
sig_regions_COM_neg <- process_hits(com_neg_q3, gr_ca)

# --- 3. Motif Discovery (CAC) ---
# Function to count absolute occurrences of 'CAC' motifs in canonical chromosomes
process_regions <- function(df, group_name) {
  gr <- GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start, end = df$end),
    mean_non_cpg_meth = df$mean_ca
  )
  
  # Restrict analysis to canonical chromosomes (1-22, X, Y, M)
  canonical_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
  gr <- gr[seqnames(gr) %in% canonical_chrs]
  
  # Extract sequences and count CAC motifs
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  n_CAC <- vcountPattern("CAC", seqs)
  
  # Format output
  df_out <- as.data.frame(gr) %>%
    mutate(
      n_CAC = n_CAC,
      group = group_name
    )
  return(df_out)
}

# --- 4. Process Target Groups and Generate Null Background ---
df_com_plus <- process_regions(sig_regions_COM_pos, "CoM +")
df_com_neg  <- process_regions(sig_regions_COM_neg, "CoM -")

# Generate size-matched background (Not Correlated) from the overarching CA dataset
set.seed(42)
df_random <- ca_data %>% sample_n(size = nrow(sig_regions_COM_pos), replace = FALSE)
df_random <- process_regions(df_random, "Not Correlated")

# --- 5. Statistical Analysis ---
# Combine all cohorts into a single dataframe for statistical modeling
df_all_box <- bind_rows(df_com_plus, df_com_neg, df_random) %>%
  mutate(group = factor(group, levels = c("CoM +", "CoM -", "Not Correlated")))

# Perform One-Way ANOVA to test for variance in CAC motifs across groups
anova_res <- aov(n_CAC ~ group, data = df_all_box)
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
      n_CAC == 0 ~ "0 motifs",
      n_CAC == 1 ~ "1 motif",
      n_CAC >= 2 ~ "≥2 motifs"
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
  filter(group == "CoM +", cat_plot %in% c("0 motifs", "≥2 motifs"))

prop_2 <- df_coMplus$prop_within_cat[df_coMplus$cat_plot == "≥2 motifs"]
prop_0 <- df_coMplus$prop_within_cat[df_coMplus$cat_plot == "0 motifs"]

ratio_2_vs_0 <- prop_2 / prop_0
cat("\nRatio (≥2 motifs / 0 motifs) for CoM +:", ratio_2_vs_0, "\n")

# --- 7. Visualization ---
# Ensure correct factor leveling for aesthetics
df_prop_plot$group <- factor(df_prop_plot$group, levels = c("CoM +", "CoM -", "Not Correlated"))

publication_plot <- ggplot(df_prop_plot, aes(x = cat_plot, y = prop_within_cat, fill = group)) +
  geom_col(alpha = 0.7)+
  scale_fill_manual(values = c("CoM +" = "#0072B2", "CoM -" = "#C22506", "Not Correlated" = "#7c807c")) +
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

# Display the plot
print(publication_plot)
