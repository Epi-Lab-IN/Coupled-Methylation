# FigS5_01_centered_bvalue_lm_panels.R
#
# Purpose
#   Generate the four-panel, mean-centered B-value linear-model figure for one
#   dataset (run this script separately for WGBS and EMseq). The only figure
#   outputs are:
#     figS5_bvalue_mean_centered_lm.png
#     figS5_bvalue_mean_centered_lm.pdf
#
# Required inputs
#   * one CG B-value window table and one CA B-value window table;
#   * both tables have chr, start, end plus one methylation column per sample.
#
# Optional sample manifest
#   If sample_manifest_file exists, it must contain: cg_column, ca_column, age,
#   sample_id. If it is NULL or missing, the legacy fallback extracts the first
#   integer from each column name. A legacy "772" age label is interpreted as
#   the second 77-year biological replicate. The resolved mapping is always
#   written to the output directory for audit.
#
# Analytical rules:
#   1. mCG strict_all_ages strata:
#      low: every sample is <33; intermediate: every sample is 33--66;
#      high: every sample is >66. A window with a missing value is unassigned.
#   2. mCA all_age_mean_cutoff_consensus strata:
#      calculate each window's all-sample mean mCA; obtain its 33rd/66th
#      percentile cutoffs across windows; then assign low/intermediate/high only
#      when every sample lies in the same cutoff-defined band.
#   3. Centering:
#      every B-value is centered by subtracting the mean across ALL joined
#      CG--CA windows in that sample, before stratum means are calculated.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(scales)
})

# ---- User paths: edit before running -----------------------------------------
dataset_label <- "WGBS"  # Set to "WGBS" or "EMseq" for the current run.
cg_bvalue_file <- "/path/to/CG_Bvalue_windows.tsv"
ca_bvalue_file <- "/path/to/CA_Bvalue_windows.tsv"
output_dir <- "/path/to/FigS5_results/WGBS"

# Optional TSV with columns: cg_column, ca_column, age, sample_id.
# Leave NULL to use legacy column-name age inference.
sample_manifest_file <- NULL

# Fixed plotting settings for the figure.
alpha <- 0.05
figure_width_in <- 12
figure_height_in <- 9.5
figure_dpi <- 400

coord_cols <- c("chr", "start", "end")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Input helpers ------------------------------------------------------------
read_bvalue_table <- function(path) {
  if (!file.exists(path)) stop("Input table not found: ", path)
  x <- readr::read_tsv(path, show_col_types = FALSE)
  missing <- setdiff(coord_cols, names(x))
  if (length(missing)) stop("Missing coordinate column(s) in ", path, ": ", paste(missing, collapse = ", "))
  if (anyDuplicated(x[, coord_cols])) stop("Coordinate rows are not unique in ", path, "; exact joining would be ambiguous.")
  x$chr <- as.character(x$chr)
  for (nm in setdiff(names(x), coord_cols)) x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
  x
}

normalise_to_percent <- function(x) {
  sample_cols <- setdiff(names(x), coord_cols)
  maximum <- max(unlist(x[sample_cols]), na.rm = TRUE)
  if (!is.finite(maximum)) stop("No finite B-values found.")
  detected_scale <- if (maximum <= 1.5) "proportion_0_1" else "percent_0_100"
  if (detected_scale == "proportion_0_1") x[sample_cols] <- lapply(x[sample_cols], function(v) v * 100)
  attr(x, "input_scale") <- detected_scale
  x
}

# Legacy fallback: extract the first integer in a sample column. A label 772 is
# a historical coding of the second 77-year replicate, not an age of 772 years.
infer_mapping <- function(x, context) {
  cols <- setdiff(names(x), coord_cols)
  raw_age <- stringr::str_match(cols, "(\\d+)")[, 2]
  age <- suppressWarnings(as.integer(raw_age))
  if (anyNA(age)) stop("Could not infer an age from every ", context, " column. Supply sample_manifest_file.")
  age[raw_age == "772"] <- 77L
  map <- tibble(input_column = cols, age = age) %>%
    group_by(age) %>% mutate(replicate = row_number(), n_replicates = n()) %>% ungroup() %>%
    mutate(sample_id = if_else(n_replicates > 1L, paste0("age", age, "_rep", replicate), paste0("age", age))) %>%
    arrange(age, replicate)
  names(map)[names(map) == "input_column"] <- paste0(context, "_column")
  map
}

resolve_mapping <- function(cg, ca, manifest_path = NULL) {
  if (!is.null(manifest_path) && file.exists(manifest_path)) {
    map <- readr::read_tsv(manifest_path, show_col_types = FALSE)
    needed <- c("cg_column", "ca_column", "age", "sample_id")
    missing <- setdiff(needed, names(map))
    if (length(missing)) stop("Sample manifest lacks: ", paste(missing, collapse = ", "))
    map <- map %>% transmute(cg_column = as.character(cg_column), ca_column = as.character(ca_column),
                             age = as.numeric(age), sample_id = as.character(sample_id)) %>% arrange(age, sample_id)
    if (anyDuplicated(map$sample_id) || anyDuplicated(map$cg_column) || anyDuplicated(map$ca_column)) stop("Manifest sample_id, cg_column, and ca_column must each be unique.")
    if (!setequal(map$cg_column, setdiff(names(cg), coord_cols))) stop("Manifest CG columns do not exactly match CG input columns.")
    if (!setequal(map$ca_column, setdiff(names(ca), coord_cols))) stop("Manifest CA columns do not exactly match CA input columns.")
    map$mapping_source <- "provided_manifest"
    return(map)
  }

  if (!is.null(manifest_path)) message("Sample manifest was not found; using legacy column-name age inference.")
  cg_map <- infer_mapping(cg, "cg")
  ca_map <- infer_mapping(ca, "ca")
  map <- dplyr::inner_join(cg_map, ca_map, by = c("age", "replicate", "n_replicates", "sample_id")) %>% arrange(age, replicate)
  if (nrow(map) != ncol(cg) - length(coord_cols) || nrow(map) != ncol(ca) - length(coord_cols)) {
    stop("CG and CA columns could not be matched by legacy age/replicate inference. Supply sample_manifest_file.")
  }
  map$mapping_source <- "legacy_column_name_inference"
  map
}

rename_and_order <- function(x, mapping, context) {
  old_cols <- mapping[[paste0(context, "_column")]]
  new_cols <- mapping$sample_id
  x <- x[, c(coord_cols, old_cols)]
  names(x) <- c(coord_cols, new_cols)
  x
}

# ---- Fixed stratum definitions ------------------------------------------------
classify_cg_strict_all_ages <- function(mat, low_threshold = 33, high_threshold = 66) {
  assignment <- apply(mat, 1, function(v) {
    if (!all(is.finite(v))) return(NA_character_)
    if (all(v < low_threshold)) return("low")
    if (all(v >= low_threshold & v <= high_threshold)) return("intermediate")
    if (all(v > high_threshold)) return("high")
    NA_character_
  })
  factor(assignment, levels = c("low", "intermediate", "high"))
}

classify_ca_all_age_mean_cutoff_consensus <- function(mat, probs = c(0.33, 0.66)) {
  window_means <- rowMeans(mat, na.rm = TRUE)
  cutoffs <- stats::quantile(window_means, probs = probs, na.rm = TRUE, names = FALSE)
  assignment <- apply(mat, 1, function(v) {
    if (!all(is.finite(v))) return(NA_character_)
    if (all(v < cutoffs[1])) return("low")
    if (all(v >= cutoffs[1] & v <= cutoffs[2])) return("intermediate")
    if (all(v > cutoffs[2])) return("high")
    NA_character_
  })
  list(strata = factor(assignment, levels = c("low", "intermediate", "high")), cut_low = cutoffs[1], cut_high = cutoffs[2])
}

stratum_counts <- function(strata, direction, total_windows) {
  lev <- c("low", "intermediate", "high")
  n <- table(factor(strata, levels = lev))
  tibble(analysis = direction, stratum = lev, n_windows = as.integer(n[lev]),
         n_unassigned = sum(is.na(strata)), total_joined_windows = total_windows)
}

# ---- Required centered mean trajectories -------------------------------------
center_by_sample_mean <- function(mat) sweep(mat, 2, colMeans(mat, na.rm = TRUE), FUN = "-")

make_trajectory_table <- function(mat, strata, mapping, analysis, measured_context, stratifying_context) {
  out <- lapply(c("high", "low"), function(group) {
    keep <- which(strata == group)
    if (!length(keep)) return(NULL)
    tibble(dataset = dataset_label, analysis = analysis, measured_context = measured_context,
           stratifying_context = stratifying_context, stratum = group,
           sample_id = mapping$sample_id, age = mapping$age,
           mean_centered_bvalue = colMeans(mat[keep, , drop = FALSE], na.rm = TRUE),
           n_windows = length(keep))
  })
  dplyr::bind_rows(out)
}

# ---- Linear-model plot helpers ------------------------------------------------
panel_theme <- function() {
  ggplot2::theme_bw(base_size = 11) + ggplot2::theme(
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = ggplot2::element_line(linewidth = 0.6),
    axis.ticks = ggplot2::element_line(linewidth = 0.8),
    axis.text = ggplot2::element_text(size = 13, color = "black"),
    axis.title = ggplot2::element_text(size = 13.5),
    plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9.5, color = "grey25"),
    panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(8, 12, 8, 14)
  )
}

fit_lm_stats <- function(df) {
  if (nrow(df) < 3 || length(unique(df$age)) < 3 || stats::sd(df$mean_centered_bvalue) == 0) {
    return(tibble(n_points = nrow(df), lm_slope = NA_real_, lm_intercept = NA_real_, lm_pvalue = NA_real_, pearson_r = NA_real_, r_squared = NA_real_))
  }
  fit <- stats::lm(mean_centered_bvalue ~ age, data = df)
  s <- summary(fit)
  pearson <- suppressWarnings(stats::cor.test(df$mean_centered_bvalue, df$age, method = "pearson"))
  tibble(n_points = nrow(df), lm_slope = unname(stats::coef(fit)[2]), lm_intercept = unname(stats::coef(fit)[1]),
         lm_pvalue = s$coefficients[2, 4], pearson_r = unname(pearson$estimate), r_squared = s$r.squared)
}

context_colours <- c(high = "#2c6fad", low = "#c0392b", not_significant = "#8d9192")

build_panel <- function(df, title, y_label) {
  stat <- fit_lm_stats(df)
  is_significant <- is.finite(stat$lm_pvalue) && stat$lm_pvalue < alpha
  colour <- if (is_significant) context_colours[[unique(df$stratum)]] else context_colours[["not_significant"]]
  subtitle <- if (is.finite(stat$lm_slope)) sprintf("LM slope = %.4g, p = %.3g | R = %.3f, R² = %.3f", stat$lm_slope, stat$lm_pvalue, stat$pearson_r, stat$r_squared) else "Linear model unavailable"
  p <- ggplot2::ggplot(df, ggplot2::aes(x = age, y = mean_centered_bvalue)) +
    ggplot2::geom_point(color = colour, size = 3, alpha = 0.85) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Age (years)", y = y_label) + panel_theme()
  if (is.finite(stat$lm_slope)) {
    p <- p + ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                                  color = "black", fill = colour, alpha = 0.20,
                                  linetype = "dashed", linewidth = 0.8)
  }
  list(plot = p, statistics = dplyr::mutate(stat, title = title, stratum = unique(df$stratum), is_significant = is_significant, colour = colour))
}

save_four_panel_figure <- function(p_tl, p_tr, p_bl, p_br, stem) {
  draw <- function() {
    grid::grid.newpage()
    layout <- grid::grid.layout(nrow = 3, ncol = 2,
      heights = grid::unit(c(1, 0.014, 1), c("null", "npc", "null")), widths = grid::unit(c(1, 1), "null"))
    grid::pushViewport(grid::viewport(layout = layout))
    print(p_tl, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_tr, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid::grid.rect(vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2), gp = grid::gpar(fill = "black", col = NA))
    print(p_bl, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
    print(p_br, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
    grid::popViewport()
  }
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  pdf_device(file.path(output_dir, paste0(stem, ".pdf")), width = figure_width_in, height = figure_height_in)
  draw(); grDevices::dev.off()
  grDevices::png(file.path(output_dir, paste0(stem, ".png")), width = figure_width_in, height = figure_height_in, units = "in", res = figure_dpi)
  draw(); grDevices::dev.off()
}

# ---- Analysis -----------------------------------------------------------------
cg_input <- normalise_to_percent(read_bvalue_table(cg_bvalue_file))
ca_input <- normalise_to_percent(read_bvalue_table(ca_bvalue_file))
mapping <- resolve_mapping(cg_input, ca_input, sample_manifest_file)
readr::write_tsv(mapping %>% mutate(dataset = dataset_label, .before = 1), file.path(output_dir, "figS5_sample_mapping_used.tsv"))

cg <- rename_and_order(cg_input, mapping, "cg")
ca <- rename_and_order(ca_input, mapping, "ca")
joined <- dplyr::inner_join(cg, ca, by = coord_cols, suffix = c("_cg", "_ca"))
if (!nrow(joined)) stop("No exact CG--CA coordinate matches were found.")

sample_ids <- mapping$sample_id
cg_mat <- as.matrix(joined[, paste0(sample_ids, "_cg")])
ca_mat <- as.matrix(joined[, paste0(sample_ids, "_ca")])
storage.mode(cg_mat) <- "numeric"; storage.mode(ca_mat) <- "numeric"

cg_strata <- classify_cg_strict_all_ages(cg_mat)
ca_result <- classify_ca_all_age_mean_cutoff_consensus(ca_mat)
ca_strata <- ca_result$strata

# Centered means used by the requested figure.
ca_centered <- center_by_sample_mean(ca_mat)
cg_centered <- center_by_sample_mean(cg_mat)
trajectories <- dplyr::bind_rows(
  make_trajectory_table(ca_centered, cg_strata, mapping, "mCA_by_mCG_strata", "mCA", "mCG"),
  make_trajectory_table(cg_centered, ca_strata, mapping, "mCG_by_mCA_strata", "mCG", "mCA")
)
if (!nrow(trajectories)) stop("No high/low windows were assigned; cannot create the requested figure.")

counts <- dplyr::bind_rows(
  stratum_counts(cg_strata, "mCA_by_mCG_strata", nrow(joined)),
  stratum_counts(ca_strata, "mCG_by_mCA_strata", nrow(joined))
)
metadata <- tibble(parameter = c("dataset", "cg_bvalue_file", "ca_bvalue_file", "sample_mapping_source", "cg_input_scale", "ca_input_scale", "join_mode", "cg_strata_rule", "cg_low_threshold", "cg_high_threshold", "ca_strata_rule", "ca_quantile_low", "ca_quantile_high", "ca_cutoff_low", "ca_cutoff_high", "joined_windows", "alpha"),
                   value = c(dataset_label, cg_bvalue_file, ca_bvalue_file, unique(mapping$mapping_source), attr(cg_input, "input_scale"), attr(ca_input, "input_scale"), "exact", "strict_all_ages", "33", "66", "all_age_mean_cutoff_consensus", "0.33", "0.66", as.character(ca_result$cut_low), as.character(ca_result$cut_high), as.character(nrow(joined)), as.character(alpha)))
readr::write_tsv(trajectories, file.path(output_dir, "figS5_centered_mean_stratum_trajectories.tsv"))
readr::write_tsv(counts, file.path(output_dir, "figS5_stratum_counts.tsv"))
readr::write_tsv(metadata, file.path(output_dir, "figS5_run_metadata.tsv"))

# ---- Four requested panels ----------------------------------------------------
get_panel_data <- function(analysis_value, stratum_value) {
  # Use argument names distinct from data-frame column names. In dplyr,
  # `analysis == analysis` or `stratum == stratum` would compare a column to
  # itself under data masking and would therefore return both strata.
  x <- trajectories %>%
    filter(.data$analysis == analysis_value, .data$stratum == stratum_value) %>%
    arrange(age, sample_id)
  if (!nrow(x)) stop("Missing required trajectory: ", analysis_value, " / ", stratum_value)
  x
}

tl <- build_panel(get_panel_data("mCA_by_mCG_strata", "high"), "mCA in high mCG context", "Mean mCA in mCG context − Mean mCA\n(B-value)")
tr <- build_panel(get_panel_data("mCA_by_mCG_strata", "low"), "mCA in low mCG context", "Mean mCA in mCG context − Mean mCA\n(B-value)")
bl <- build_panel(get_panel_data("mCG_by_mCA_strata", "high"), "mCG in high mCA context", "Mean mCG in mCA context − Mean mCG\n(B-value)")
br <- build_panel(get_panel_data("mCG_by_mCA_strata", "low"), "mCG in low mCA context", "Mean mCG in mCA context − Mean mCG\n(B-value)")

panel_stats <- dplyr::bind_rows(
  dplyr::mutate(tl$statistics, position = "top_left", analysis = "mCA_by_mCG_strata"),
  dplyr::mutate(tr$statistics, position = "top_right", analysis = "mCA_by_mCG_strata"),
  dplyr::mutate(bl$statistics, position = "bottom_left", analysis = "mCG_by_mCA_strata"),
  dplyr::mutate(br$statistics, position = "bottom_right", analysis = "mCG_by_mCA_strata")
)
readr::write_tsv(panel_stats, file.path(output_dir, "figS5_bvalue_mean_centered_lm_panel_statistics.tsv"))

figure_stem <- "figS5_panel_bvalue_mean_centered_lm"
save_four_panel_figure(tl$plot, tr$plot, bl$plot, br$plot, figure_stem)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
message("Complete: ", file.path(output_dir, paste0(figure_stem, ".png")))
