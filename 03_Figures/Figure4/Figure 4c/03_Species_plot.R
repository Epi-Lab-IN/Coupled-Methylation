# ===============================================================
# Cross-species comparison of CoM+ CA methylation windows
# ===============================================================

# -------------------------
# 0. Load required libraries
# -------------------------
library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(stringr)

in_dir <- "/path/to/Fig4c/02_species_CoM_correlation/"

# -------------------------
# 1. Read data for each species and add species labels
# -------------------------
platypus_df <- read_delim(
  paste0(in_dir,"cowindows_platypus.tsv"),
  delim = "\t"
) %>%
  mutate(species = "Platypus",  chr = as.character(chr))

chicken_df <- read_delim(
  paste0(in_dir,"cowindows_chicken.tsv"),
  delim = "\t"
) %>%
  mutate(species = "Chicken", chr = as.character(chr))

frog_df <- read_delim(
  paste0(in_dir,"cowindows_frog.tsv"),
  delim = "\t"
) %>%
  mutate(species = "Frog", chr = as.character(chr))

zebrafish_df <- read_delim(
  paste0(in_dir,"cowindows_zebrafish.tsv"),
  delim = "\t"
) %>%
  mutate(species = "Zebrafish", chr = as.character(chr))

shark_df <- read_delim(
  paste0(in_dir,"cowindows_shark.tsv"),
  delim = "\t"
) %>%
  mutate(species = "Elephant shark", chr = as.character(chr))

# -------------------------
# 2. Combine all species into one dataframe
# -------------------------
df_all_species <- bind_rows(platypus_df, chicken_df, frog_df, zebrafish_df, shark_df)

# -------------------------
# 3. Filter high-methylation CA windows (top 25%) per species
# -------------------------
df_filtered <- df_all_species %>%
  group_by(species) %>%
  filter(CA_bvalues > quantile(CA_bvalues, 0.66, na.rm = TRUE)) %>%
  ungroup()

# -------------------------
# 4. Create bins for mCG β-values (0–100%, step = 10)
# -------------------------
df_filtered <- df_filtered %>%
  mutate(CG_bin = cut(CG_bvalues, breaks = seq(0, 100, by = 10), include.lowest = TRUE))

# -------------------------
# 5. Count number of high-mCA events per CG bin and species
# -------------------------
df_counts <- df_filtered %>%
  group_by(species, CG_bin) %>%
  summarise(n_events = n(), .groups = "drop") %>%
  mutate(
    bin_nums = str_extract_all(CG_bin, "\\d+\\.?\\d*"),
    bin_lower = as.numeric(sapply(bin_nums, `[`, 1)),
    bin_upper = as.numeric(sapply(bin_nums, `[`, 2)),
    CG_bin_center = (bin_lower + bin_upper) / 2
  ) %>%
  dplyr::select(-bin_nums, -bin_lower, -bin_upper)

# -------------------------
# 6. Normalize event counts by genome size (in Gb)
# -------------------------
genome_sizes <- data.frame(
  species = c("Platypus", "Chicken", "Frog", "Zebrafish", "Elephant shark"),
  genome_size = c(1.8, 1.0, 1.4, 1.3, 0.91)  # approximate genome sizes (Gb)
)

df_counts <- df_counts %>%
  left_join(genome_sizes, by = "species") %>%
  mutate(n_events_norm = n_events / genome_size)

# -------------------------
# 7. Prepare data for plotting
# -------------------------
# Normalize further to match scale (×10e5)
df_counts$n_events_norm <- df_counts$n_events_norm / 100000

# Okabe-Ito color palette (colorblind-friendly)
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7"
)

# -------------------------
# 8. Plot species-specific trends in CA vs CG methylation
# -------------------------
ggplot(df_counts, aes(x = CG_bin_center, y = n_events_norm, color = species)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.7, linewidth = 2) +
  scale_x_continuous(
    breaks = seq(0, 100, 25),
    limits = c(0, 100)
  ) +
  scale_color_manual(values = okabe_ito) +
  labs(
    x = "mCG (β-value)",
    y = "Number of mCA windows (×10⁵)",
    color = "Species"
  ) +
  theme_classic(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 18),
    axis.ticks = element_line(size = 1),
    axis.line = element_line(size = 1),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  ) 
