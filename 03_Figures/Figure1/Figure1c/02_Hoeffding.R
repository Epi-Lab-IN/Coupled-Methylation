# ================================
# Load required libraries
# ================================
library(readr)
library(GenomicRanges)
library(Hmisc)
library(ggplot2)
library(patchwork)
library(scales)

# ================================
# Helper Functions
# ================================

# Read a tab-delimited file and convert it to a GRanges object
read_to_granges <- function(path) {
  df <- read_delim(path, delim = "\t")
  makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
}

# Merge overlapping regions between CA and CG contexts
merge_overlaps <- function(cg_GR_input, ca_GR_input) {
  overlaps <- findOverlaps(ca_GR_input, cg_GR_input)
  
  cg_df <- as.data.frame(cg_GR_input)
  ca_df <- as.data.frame(ca_GR_input)
  
  merged <- cbind(
    cg_df[subjectHits(overlaps), ],
    ca_df[queryHits(overlaps), ]
  )
  
  merged$region_id <- paste(
    gsub("chr", "", merged[, 1]), 
    merged[, 2],                  
    merged[, 3],                  
    sep = "-"
  )
  
  row.names(merged) <- merged$region_id
  
  return(merged)
}

# Compute the Hoeffding test on the column means of a subset
compute_hoeffding <- function(subset_matrix, ages) {
  mean_values <- colMeans(subset_matrix, na.rm = TRUE)
  hoeffd(mean_values, ages)
}

# Extract age values from column names (e.g., "Sample_25y" -> 25)
extract_ages <- function(col_names) {
  as.numeric(sub(".*?(\\d+)y$", "\\1", col_names))
}

# Extract relevant stats (D and P values) from a Hoeffding test object
get_hoeff_stats <- function(hoeff_obj) {
  d_val <- hoeff_obj$D[1, 2]
  p_val <- hoeff_obj$P[1, 2]
  list(D = d_val, P = p_val)
}

# Perform Hoeffding tests and build the data frames used for plotting
analyze_correlation <- function(data, cg_cols, ca_cols, threshold_type = "cg") {
  col_names_active <- colnames(data)[ca_cols]
  ages <- extract_ages(col_names_active)
  
  res <- list(stats = list(), plot_data = list())
  
  if (threshold_type == "cg") {
    # Threshold based on CG context
    high_set <- data[apply(data[, cg_cols], 1, function(x) all(x > 66, na.rm = TRUE)), ]
    low_set  <- data[apply(data[, cg_cols], 1, function(x) all(x < 33, na.rm = TRUE)), ]
    
    res$stats$high <- compute_hoeffding(high_set[, ca_cols], ages)
    res$stats$low  <- compute_hoeffding(low_set[, ca_cols], ages)
    
    means_high <- colMeans(high_set[, ca_cols], na.rm = TRUE)
    means_low  <- colMeans(low_set[, ca_cols], na.rm = TRUE)
    
    res$plot_data$high <- data.frame(Age = ages, MeanMeth = means_high, Group = "High mCG", Context = "mCA")
    res$plot_data$low  <- data.frame(Age = ages, MeanMeth = means_low,  Group = "Low mCG",  Context = "mCA")
    
  } else {
    # Threshold based on CA context
    mean_ca_values <- rowMeans(data[, ca_cols], na.rm = TRUE)
    bot <- quantile(mean_ca_values, 0.33, na.rm = TRUE)
    top <- quantile(mean_ca_values, 0.66, na.rm = TRUE)
    
    high_set <- data[apply(data[, ca_cols], 1, function(x) all(x > top, na.rm = TRUE)), ]
    low_set  <- data[apply(data[, ca_cols], 1, function(x) all(x < bot, na.rm = TRUE)), ]
    
    res$stats$high <- compute_hoeffding(high_set[, cg_cols], ages)
    res$stats$low  <- compute_hoeffding(low_set[, cg_cols], ages)
    
    means_high <- colMeans(high_set[, cg_cols], na.rm = TRUE)
    means_low  <- colMeans(low_set[, cg_cols], na.rm = TRUE)
    
    ages_cg <- extract_ages(colnames(data)[cg_cols])
    
    res$plot_data$high <- data.frame(Age = ages_cg, MeanMeth = means_high, Group = "High mCA", Context = "mCG")
    res$plot_data$low  <- data.frame(Age = ages_cg, MeanMeth = means_low,  Group = "Low mCA",  Context = "mCG")
  }
  
  return(res)
}

# Build an individual ggplot, optionally with fixed y-axis limits and breaks
create_individual_plot <- function(df_plot, title_str, y_lab, color_point, y_limits = NULL, y_breaks = NULL) {
  p <- ggplot(df_plot, aes(x = Age, y = MeanMeth)) +
    geom_smooth(method = "loess", span = 2, se = TRUE, color = "black", fill = color_point, alpha = 0.2, linetype = "dashed") +
    geom_point(color = color_point, size = 3, alpha = 0.8) +
    theme_bw(base_size = 10) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_line(linewidth = 1),  
      axis.ticks = element_line(linewidth = 1), 
      axis.text = element_text(size = 10, color = "black"),   
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold") 
    ) +
    labs(
      title = title_str,
      subtitle = df_plot$Group[1],
      x = "Age (years)",
      y = y_lab
    )
  
  # coord_cartesian only zooms visually; it does not drop data/CI outside the range
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  if (!is.null(y_breaks)) {
    p <- p + scale_y_continuous(breaks = y_breaks, labels = number_format(accuracy = 0.1))
  }
  
  p
}

# Get the y-axis range ggplot will actually draw (includes the geom_smooth CI band)
get_rendered_y_range <- function(df_plot) {
  p_tmp <- ggplot(df_plot, aes(x = Age, y = MeanMeth)) +
    geom_smooth(method = "loess", span = 2, se = TRUE)
  built <- ggplot_build(p_tmp)
  built$layout$panel_params[[1]]$y.range
}

# Match two y-axis ranges to the same width, extending the shorter one
match_y_ranges <- function(range_a, range_b, extend_side = "center") {
  width_a <- diff(range_a)
  width_b <- diff(range_b)
  target_width <- max(width_a, width_b)
  
  adjust <- function(rng) {
    needed <- target_width - diff(rng)
    if (needed <= 0) return(rng)
    if (extend_side == "top") {
      c(rng[1], rng[2] + needed)
    } else if (extend_side == "bottom") {
      c(rng[1] - needed, rng[2])
    } else { # "center"
      c(rng[1] - needed / 2, rng[2] + needed / 2)
    }
  }
  
  list(a = adjust(range_a), b = adjust(range_b))
}

# Generate identical breaks (same spacing) for two ranges that already share the same width
get_matched_breaks <- function(range_a, range_b, n_breaks = 3) {
  width <- diff(range_a)  # range_a and range_b already share the same width
  list(
    a = seq(range_a[1], range_a[1] + width, length.out = n_breaks),
    b = seq(range_b[1], range_b[1] + width, length.out = n_breaks)
  )
}

# ================================
# Data index configuration
# ================================
N_INFO_COLS <- 5
N_SAMPLES <- 7
idx_cg_data <- (N_INFO_COLS + 1) : (N_INFO_COLS + N_SAMPLES)
idx_ca_data <- (N_INFO_COLS + N_SAMPLES + N_INFO_COLS + 1) : (N_INFO_COLS + N_SAMPLES + N_INFO_COLS + N_SAMPLES)

# ================================
# mCA Context Analysis
# ================================
cg_path_bval <- "/home/epi-lab/Escritorio/Resultados/CHH/EMseq_cg_comunes_bvalues_cowindows.tsv"
ca_path_mval <- "/home/epi-lab/Escritorio/Resultados/CHH/EMseq_ca_comunes_Mvalues_cowindows.tsv"

cg_b <- read_to_granges(cg_path_bval)
ca_m <- read_to_granges(ca_path_mval)

cg_ca_combined <- merge_overlaps(cg_b, ca_m)

cg_results_full <- analyze_correlation(
  data = cg_ca_combined,
  cg_cols = idx_cg_data,
  ca_cols = idx_ca_data,
  threshold_type = "cg"
)

print("mCA statistics:")
print(cg_results_full$stats)

# ================================
# mCG Context Analysis
# ================================
cg_path_mval_2 <- "/home/epi-lab/Escritorio/Resultados/CHH/EMseq_cg_comunes_Mvalues_cowindows.tsv"
ca_path_bval_2 <- "/home/epi-lab/Escritorio/Resultados/CHH/EMseq_ca_comunes_bvalues_cowindows.tsv"

cg_m_2 <- read_to_granges(cg_path_mval_2)
ca_b_2 <- read_to_granges(ca_path_bval_2)

ca_cg_combined <- merge_overlaps(cg_m_2, ca_b_2)

ca_results_full <- analyze_correlation(
  data = ca_cg_combined,
  cg_cols = idx_cg_data,
  ca_cols = idx_ca_data,
  threshold_type = "ca"
)

print("mCG statistics:")
print(ca_results_full$stats)

# ================================
# 4-Panel Plot Generation (fixed span = 2)
# ================================
y_lab_mca <- "Mean mCA (M-values)"
y_lab_mcg <- "Mean mCG (M-values)"
loess_span <- 2

# Match y-axis ranges (including CI) within each pair of plots
range_high_mCG <- get_rendered_y_range(cg_results_full$plot_data$high)
range_low_mCG  <- get_rendered_y_range(cg_results_full$plot_data$low)
matched_mCA <- match_y_ranges(range_high_mCG, range_low_mCG, extend_side = "center")
breaks_mCA  <- get_matched_breaks(matched_mCA$a, matched_mCA$b)

range_high_mCA <- get_rendered_y_range(ca_results_full$plot_data$high)
range_low_mCA  <- get_rendered_y_range(ca_results_full$plot_data$low)
matched_mCG <- match_y_ranges(range_high_mCA, range_low_mCA, extend_side = "center")
breaks_mCG  <- get_matched_breaks(matched_mCG$a, matched_mCG$b)

# Build the four individual plots
p1 <- create_individual_plot(cg_results_full$plot_data$high, "mCA vs Age", y_lab_mca, "#0072b2", y_limits = matched_mCA$a, y_breaks = breaks_mCA$a)
p2 <- create_individual_plot(cg_results_full$plot_data$low,  "mCA vs Age", y_lab_mca, "#7c807c", y_limits = matched_mCA$b, y_breaks = breaks_mCA$b)

p3 <- create_individual_plot(ca_results_full$plot_data$high, "mCG vs Age", y_lab_mcg, "#7c807c",  y_limits = matched_mCG$a, y_breaks = breaks_mCG$a)
p4 <- create_individual_plot(ca_results_full$plot_data$low,  "mCG vs Age", y_lab_mcg, "#C22506",  y_limits = matched_mCG$b, y_breaks = breaks_mCG$b)

# Combine into a single 2x2 panel
final_panel <- (p1 | p2) / (p3 | p4)

print(final_panel)