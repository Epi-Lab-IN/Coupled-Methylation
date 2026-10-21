# ================================================================
# LIBRARIES & SETUP
# ================================================================
suppressPackageStartupMessages({
  library(readr)
  library(data.table)  # For fast reading (fread)
  library(dplyr)       # For data manipulation
  library(ggplot2)     # For plotting
  library(broom)       # For tidy model summaries
  library(matrixStats) # For ultra-fast row medians
})

# ================================================================
# 1. LOAD SIGNIFICANT GENE SETS
# ================================================================
genes_pos <- read_delim("~/Escritorio/Resultados/CHH/genes/genes_COM_pos_aggregated.txt", 
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

genes_neg <- read_delim("~/Escritorio/Resultados/CHH/genes/genes_COM_neg_aggregated.txt", 
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# ================================================================
# 2. LOAD & FORMAT RNA-SEQ / MICROARRAY DATASETS
# ================================================================
data <- fread("/home/epi-lab/Escritorio/RNA_GSE/3_RNA_datasets_updated.txt", header = TRUE)

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
    
    # Protect against division by zero in entirely empty columns
    if (total_valid > 0 && (valid_numeric / total_valid) > 0.95) {
      return(x_num)
    } else {
      return(x)
    }
  }
  return(x)
})

# ================================================================
# 3. PLOTTING FUNCTION
# ================================================================
# Note: Added 'plot_title' argument to differentiate CoM+ and CoM-
plot_median_expression_noCollapse <- function(genes, controles_df, plot_color, plot_title) {
  
  # --- Filter for valid genes present in the dataset ---
  genes_rna <- intersect(colnames(controles_df), genes)
  genes_rna <- genes_rna[sapply(controles_df[genes_rna], is.numeric)]
  
  if(length(genes_rna) == 0){
    stop("No valid numeric genes found to process.")
  }
  
  # =========================
  # A. TARGET GENES (REAL)
  # =========================
  rna_sel <- controles_df %>%
    dplyr::select(Age, all_of(genes_rna))
  
  # matrixStats::rowMedians is exponentially faster than apply(..., 1, median)
  rna_sel_matrix <- as.matrix(rna_sel %>% dplyr::select(-Age))
  rna_sel$median_expr <- rowMedians(rna_sel_matrix, na.rm = TRUE)
  
  modelo_real <- lm(median_expr ~ Age, data = rna_sel)
  r2_real <- summary(modelo_real)$r.squared
  p_real  <- broom::tidy(modelo_real)$p.value[2]
  
  # Spearman correlation (non-parametric)
  spearman_real <- cor.test(rna_sel$Age, rna_sel$median_expr, method = "spearman")
  rho_real <- spearman_real$estimate
  p_spearman_real <- spearman_real$p.value
  
  # =========================
  # B. RANDOM GENES (NO OVERLAP)
  # =========================
  set.seed(1234) # Ensure reproducibility for the random background
  n_genes <- length(genes_rna)
  
  # Define the universe of numeric genes
  universe_genes <- names(controles_df)[sapply(controles_df, is.numeric)]
  
  # Remove Age and the target genes to prevent overlap
  universe_genes <- setdiff(universe_genes, c("Age", genes_rna))
  
  if(length(universe_genes) < n_genes){
    stop("Not enough genes to construct a random non-overlapping background.")
  }
  
  genes_random <- sample(universe_genes, n_genes)
  
  rna_random <- controles_df %>%
    dplyr::select(Age, all_of(genes_random))
  
  rna_random_matrix <- as.matrix(rna_random %>% dplyr::select(-Age))
  rna_random$median_expr <- rowMedians(rna_random_matrix, na.rm = TRUE)
  
  modelo_random <- lm(median_expr ~ Age, data = rna_random)
  r2_random <- summary(modelo_random)$r.squared
  p_random  <- broom::tidy(modelo_random)$p.value[2]
  
  # Spearman correlation for random background
  spearman_random <- cor.test(rna_random$Age, rna_random$median_expr, method = "spearman")
  rho_random <- spearman_random$estimate
  p_spearman_random <- spearman_random$p.value
  
  # --- Print Stats to Console ---
  cat("\n==============================\n")
  cat("Target Genes -> R² =", round(r2_real, 3), " | P (lm) =", signif(p_real, 3), "\n")
  cat("Target Genes -> Spearman rho =", round(rho_real, 3), " | P (spearman) =", signif(p_spearman_real, 3), "\n")
  cat("Random Genes -> R² =", round(r2_random, 3), "| P (lm) =", signif(p_random, 3), "\n")
  cat("Random Genes -> Spearman rho =", round(rho_random, 3), " | P (spearman) =", signif(p_spearman_random, 3), "\n")
  cat("==============================\n")
  
  # --- Visualization ---
  g <- ggplot() +
    # Background (Random)
    geom_point(data = rna_random, aes(x = Age, y = median_expr),
               color = "grey70", alpha = 0.2, size = 2) +
    geom_smooth(data = rna_random, aes(x = Age, y = median_expr),
                method = "lm", se = TRUE, color = "grey70", fill = "grey90") +
    # Foreground (Real)
    geom_point(data = rna_sel, aes(x = Age, y = median_expr),
               color = plot_color, alpha = 0.2, size = 2) +
    geom_smooth(data = rna_sel, aes(x = Age, y = median_expr),
                method = "lm", se = TRUE, color = plot_color, fill = plot_color, alpha = 0.15) +
    
    labs(x = "Age", y = "RNA Z-score", title = plot_title) +
    coord_cartesian(ylim = c(-0.5, 0.5)) + 
    scale_x_continuous(breaks = seq(40, 100, by = 20)) +
    
    theme_bw(base_size = 14) +
    theme(
      
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_line(linewidth = 1),  
      axis.ticks = element_line(linewidth = 1), 
      axis.text = element_text(size = 18, color = "black"),   
      axis.title = element_text(size = 19),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 19),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold") 
    )
  
  print(g)
}

# ================================================================
# 4. EXECUTE PLOTS
# ================================================================
# Plot CoM+
plot_median_expression_noCollapse(
  genes = genes_pos$gene_symbol,
  controles_df = controles,
  plot_color = "#4a0480",
  plot_title = "CoM+ genes (GSE36192+GSE33000+GSE15222)"
)

# Plot CoM-
plot_median_expression_noCollapse(
  genes = genes_neg$gene_symbol, 
  controles_df = controles, 
  plot_color = "#ab8f02",
  plot_title = "CoM- genes (GSE36192+GSE33000+GSE15222)"
)