# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
## Plot the correlation between CG and CA methylation using 2D density

# Define a custom grayscale palette
greys <- colorRampPalette(c("#ffffff", "#3d3d3d", "#1f1e1e", "#0a0a0a"))(100)

# Create a filled 2D density plot of CG vs CA methylation
ggplot(cowindows, aes(x = CG_methylation, y = CA_methylation)) +
  geom_density_2d_filled(bins = 100) +  # Filled density plot with many bins
  scale_fill_manual(values = greys) +   # Apply custom grayscale
  theme_classic() +                     # Use clean theme
  theme(legend.position = "none") +     # Hide legend
  ylim(0, 40) +                          # Limit y-axis
  labs(
    x = "% mCG/CG",
    y = "% mCA/CA"
  ) +
  coord_fixed(ratio = 2.5)              # Set aspect ratio for better visual proportion
## Create methylation heatmaps with ComplexHeatmap


# Set seed for reproducibility of sampling
set.seed(123)

# Randomly sample 5000 windows from each methylation category (unmethylated, hemimethylated, methylated)
sub_datanomet <- UMwindows[sample(1:nrow(UMwindows), size = 5000, replace = FALSE), ]
sub_datasemimet <- HMwindows[sample(1:nrow(HMwindows), size = 5000, replace = FALSE), ]
sub_datamet <- Mwindows[sample(1:nrow(Mwindows), size = 5000, replace = FALSE), ]

# Combine all subsets into a single dataset
windows <- rbind(
  sub_datamet,
  sub_datasemimet,
  sub_datanomet
)

# Extract CG and CA methylation values into a data frame
data <- data.frame(
  cg = windows$CG_methylation,
  ca = windows$CA_methylation
)

# Sort rows by CG methylation in descending order for visual alignment
data <- data[order(data$cg, decreasing = TRUE), ]

# Calculate CA methylation percentiles to define color breaks
percentiles_ca <- quantile(
  data$ca,
  probs = c(0, 0.3, 0.6, 0.8)
)

# Create a color scale for CG methylation (low to high)
col_cg <- colorRamp2(
  c(min(data$cg), max(data$cg)),  
  c("#fcfcfa", "#c22506")
)

# Create a color scale for CA methylation using custom percentiles
col_ca <- colorRamp2(
  c(
    percentiles_ca[1],
    percentiles_ca[2],
    percentiles_ca[3],
    percentiles_ca[4]
  ),  
  c("white", "white", "#bad3e0", "#0072b2")
)

# Generate heatmap for CG methylation
bar_cg <- Heatmap(
  data$cg,
  name = "CG",
  col = col_cg,
  width = unit(10, "mm"),
  cluster_rows = FALSE,         # Disable clustering
  show_row_names = FALSE,       # Hide row names for clarity
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 12),
    title_gp = gpar(fontsize = 12)
  )
)

# Generate heatmap for CA methylation
bar_ca <- Heatmap(
  data$ca,
  name = "CA",
  col = col_ca,
  width = unit(10, "mm"),
  cluster_rows = FALSE,
  show_row_names = FALSE,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 12),
    title_gp = gpar(fontsize = 12)
  )
)

# Combine both heatmaps side by side
heatmap_list <- bar_cg + bar_ca
