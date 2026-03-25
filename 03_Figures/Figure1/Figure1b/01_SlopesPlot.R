library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
# --- Step 0: Read data file  ---
df <-read_delim("~/Escritorio/Paper_scripts/Figure1/Figura1b/Emseqslopes.txt", 
                                delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)


# --- Step 1: Linear expansion of each sample from x = -15 to 15 ---
expanded_df <- df %>%
  group_by(Sample) %>%
  do({
    x_vals <- seq(-15, 15)
    data.frame(
      Sample = .$Sample[1],
      Slope = .$Slope[1],
      Intercept = .$Intercept[1],
      x = x_vals,
      y = .$Slope[1] * x_vals + .$Intercept[1]
    )
  }) %>%
  ungroup()

# --- Step 2: Normalize y values so that all start from 0 at x = -15 ---
normalized_df <- expanded_df %>%
  group_by(Sample) %>%
  mutate(
    y_at_minus15 = y[x == -15][1],
    y_norm = y + abs(y_at_minus15)
  ) %>%
  ungroup()

# --- Step 3: Extract normalized y values for the "2y" sample ---
df_2y <- normalized_df %>%
  filter(Sample == "35y") %>%
  dplyr::select(x, y_2y = y_norm)

# --- Step 4: Relativize all samples to "2y" sample at each x ---
rel_df <- normalized_df %>%
  left_join(df_2y, by = "x") %>%
  mutate(y_rel = y_norm - y_2y)

# --- Step 5: Extract numeric Age from Sample name (e.g., "77y" -> 77) ---
rel_df <- rel_df %>%
  mutate(Age = as.numeric(gsub("y.*", "", Sample)))

# --- Step 6: Prepare data for plotting ---
plot_df <- rel_df %>%
  dplyr::select(Position = x, Sample, Age, Fraction = y_rel)

# --- Step 7: Create color gradient (darker = older) ---
# Ensure Age is treated as factor with correct order
plot_df <- plot_df %>%
  arrange(Age) %>%
  mutate(Age = factor(Age, levels = sort(unique(Age), decreasing = FALSE)))

# --- Step 9: Create blue color palette and order Age factor descending ---

# Create blue color palette from dark (older) to light (younger)
blue_palette <- colorRampPalette(c("#035887", "#c6dbef"))(length(unique(plot_df$Age)))

# Arrange Age ascending and set factor levels in ascending order (youngest to oldest)
plot_df <- plot_df %>%
  arrange(Age) %>%
  mutate(Age = factor(Age, levels = sort(unique(Age), decreasing = TRUE)))  # ascending order


# --- Step 10: Generate the plot ---
ggplot(plot_df, aes(x = Position, y = Fraction, group = Sample, color = Age)) +
  geom_line(size = 2) +  # Línea más gruesa
  
  # Y-axis configuration
  scale_y_continuous(
    limits = c(-0.1, 1.25),
    breaks = seq(0, 1.25, by = 0.3125),
    expand = c(0, 0.001),labels = number_format(accuracy = 0.1)  # pequeño padding
  ) +
  
  # X-axis configuration
  scale_x_continuous(breaks = seq(-15, 15, by =7.5)) +
  
  # Apply manual color palette
  scale_color_manual(values = blue_palette) +
  
  # Set theme and label elements
  theme_classic(base_size = 18) + # Texto más grande
  labs(
    x = "mCG (M-value)",
    y = "Rel. ΔmCA (M-value)",
    color = "Age"
  ) +
  theme(
    legend.position = "right",
    axis.line = element_line(size = 1),  # Grosor de los ejes
    axis.ticks = element_line(size = 1), # Grosor de las marcas
    axis.text = element_text(size = 18, color = "black"),   # Tamaño de números de ejes
    axis.title = element_text(size = 19),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 19),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  # Fixed aspect ratio
  coord_fixed(ratio = 20)
