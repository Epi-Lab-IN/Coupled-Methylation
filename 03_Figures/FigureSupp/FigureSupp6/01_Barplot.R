# ===============================================================
# Summary of CoM+ and CoM– genomic windows and visualization
# ===============================================================

# -------------------------
# 0. Load required libraries
# -------------------------
library(readr)
library(dplyr)
library(ggplot2)

# -------------------------
# 1. Read input data
# -------------------------
# CoM+ windows (positively correlated with age)
com_pos_df <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/Corr_with_age/COM_pos_filtered",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE
)

# Create unique genomic location identifier
com_pos_df$loc <- paste(gsub("chr", "", com_pos_df$chr), com_pos_df$start, com_pos_df$end, sep = "_")

# CoM– windows (negatively correlated with age)
com_neg_df <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/Corr_with_age/COM_neg_filtered",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE
)

# Create unique genomic location identifier
com_neg_df$loc <- paste(gsub("chr", "", com_neg_df$chr), com_neg_df$start, com_neg_df$end, sep = "_")

# -------------------------
# 2. Build summary table
# -------------------------
# Count unique genomic windows per group
df_counts <- tibble(
  group = c("CoM +", "CoM -"),
  count = c(length(unique(com_pos_df$loc)),
            length(unique(com_neg_df$loc)))
) %>%
  mutate(
    group = factor(group, levels = c("CoM +", "CoM -")),
    count = count / 10000  # scale for plotting
  )

# -------------------------
# 3. Plot bar chart
# -------------------------
ggplot(df_counts, aes(x = group, y = count, fill = group)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("CoM +" = "#0072B2",
                               "CoM -" = "#C22506")) +
  labs(
    x = NULL,
    y = "Number of windows (×10⁴)"
  ) +
  theme_linedraw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
