# ================================================================
# LIBRARIES & SETUP
# ================================================================
suppressPackageStartupMessages({
  library(data.table)  # For ultra-fast data loading
  library(dplyr)       # For data manipulation
  library(tidyr)       # For pivot_longer (wide to long formatting)
  library(ggplot2)     # For plotting
})

# ================================================================
# 1. LOAD & FORMAT RNA-SEQ / MICROARRAY DATASETS
# ================================================================
data <- fread("/home/epi-lab/Escritorio/RNA_GSE/3_datasets_updated.txt", header = TRUE)

# Transpose matrix: set first column as names, then remove it
controles <- t(data[, -1, with = FALSE])
colnames(controles) <- data$V1
controles <- as.data.frame(controles, stringsAsFactors = FALSE)

# Clean strings and convert columns to numeric dynamically
controles[] <- lapply(controles, function(x) {
  if (is.character(x)) {
    x <- trimws(x) # Remove whitespace
    x[x %in% c("NA", "NaN", "", ".", "Inf", "-Inf")] <- NA
    
    x_num <- suppressWarnings(as.numeric(x))
    
    # If >= 95% of non-NA values are numeric, convert the whole column
    valid_numeric <- sum(!is.na(x_num))
    total_valid   <- sum(!is.na(x))
    
    if (total_valid > 0 && (valid_numeric / total_valid) > 0.95) {
      return(x_num)
    } else {
      return(x)
    }
  }
  return(x)
})

# ================================================================
# 2. HEATMAP FUNCTION
# ================================================================
plot_rna_heatmap <- function(target_genes, data_df, high_color) {
  
  # --- Safety Check: Ensure target genes exist in the dataframe ---
  valid_genes <- intersect(target_genes, colnames(data_df))
  
  if (length(valid_genes) < length(target_genes)) {
    missing <- setdiff(target_genes, valid_genes)
    warning("The following genes were not found in the dataset and will be excluded: ", 
            paste(missing, collapse = ", "))
  }
  
  if (length(valid_genes) == 0) {
    stop("None of the target genes were found in the dataset.")
  }
  
  # --- Data Processing ---
  df_plot <- data_df %>%
    # A. Select only Age and the valid target genes
    dplyr::select(Age, all_of(valid_genes)) %>%
    
    # B. Reshape from WIDE (Columns) to LONG (Rows) format
    pivot_longer(
      cols = -Age,             # All columns except Age are genes
      names_to = "Gene",       # New column for gene names
      values_to = "Z_Score"    # New column for expression values
    ) %>%
    
    # C. Create age groups (10-year bins)
    # Breaks at 30, 40... generate centers at 35, 45...
    mutate(
      Age_Bin = cut(Age, 
                    breaks = c(30, 40, 50, 60, 70, 80, 90, 110), 
                    labels = c("35 years", "45 years", "55 years", "65 years", 
                               "75 years", "85 years", "103 years"),
                    include.lowest = TRUE)
    ) %>%
    
    # D. Calculate the median Z-score for each Gene + Age Bin combination
    group_by(Gene, Age_Bin) %>%
    summarise(Median_Z = median(Z_Score, na.rm = TRUE), .groups = "drop")
  
  # --- Ensure correct Y-axis gene order (matches input vector) ---
  df_plot$Gene <- factor(df_plot$Gene, levels = rev(target_genes))
  
  # --- Visualization ---
  p <- ggplot(df_plot, aes(x = Age_Bin, y = Gene, fill = Median_Z)) +
    # Use tiles with no borders to stitch smoothly to adjacent plots
    geom_tile(width = 1, height = 1) + 
    
    scale_fill_gradient(
      low = "#F0F0F0",     # Light grey background
      high = high_color,   # Dynamic color passed to function
      name = "RNA Z-score"
    ) +
    
    labs(x = NULL, y = NULL) + 
    
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_blank(), # Hidden so it can be stitched alongside a methylation plot
      panel.grid = element_blank(),
      legend.title = element_text(size = 14),
      plot.margin = margin(l = 10, r = 0, t = 0, b = 0)
    )
  
  return(p)
}

# ================================================================
# 3. EXECUTE PLOT
# ================================================================
# Define the exact order of genes (must match your methylation plot)
genes_to_plot <- c("APP", "CELF2", "IDE", "PRDM2") 

# Generate the plot
# You can change high_color to "#4a0480" for CoM+ or "#ab8f02" for CoM-
 plot_rna_heatmap(
  target_genes = genes_to_plot, 
  data_df = controles, 
  high_color = "#ab8f02" 
)