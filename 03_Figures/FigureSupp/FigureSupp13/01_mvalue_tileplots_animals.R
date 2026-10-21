#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## Methylation Analysis: M-value Tile Plots (animals)
##
## Description:
##   Generates 2D density tile plots of genomic methylation
##   M-values (mCG vs mCA) along with count-density contours, linear regression
##   models, and beta-to-M reference coordinate overlay markers.
##
## Environment:
##  Executed in R version 4.6.1.
##  Ubuntu 20.04.6 LTS
##
## Code Structure:
##   1. USER CONFIGURATION (Top Section):
##      All user-tunable parameters, including input/output directories, sample
##      label overrides, beta-to-M coverage reference settings, grid binning,
##      colour palette scaling, contour settings, regression options, and export
##      dimensions.
##   2. EXECUTION PIPELINE & HELPER FUNCTIONS (Bottom Section):
##      Modular functions for data ingestion, reference coordinate transformation,
##      density binning, 2D Gaussian contour smoothing, linear regression modeling,
##      ggplot2 visualization, and dual-format (PNG + Cairo PDF) figure export.
## -----------------------------------------------------------------------------

library(data.table)
library(ggplot2)

## =============================================================================
## USER CONFIGURATION ----
## =============================================================================

## ---- 1. File Paths & Dataset Settings ---------------------------------------
## Input directory containing per-sample mCG + mCA window files ('cowindows_*.tsv')
input_dir <- file.path("path/to/data/directory", "windows_animals")

## Output directory for exported figures (PNG raster and Cairo PDF vector) and tables
out_dir <- file.path("path/to/output/directory", "mvalue_tileplots_animals")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Human-readable label overrides for sample plot titles
sample_label_overrides <- c(
  "cowindows_chicken.tsv"   = "Chicken",
  "cowindows_frog.tsv"      = "Frog",
  "cowindows_platypus.tsv"  = "Platypus",
  "cowindows_shark.tsv"     = "Elephant Shark",
  "cowindows_zebrafish.tsv" = "Zebrafish"
)

## ---- 2. Beta-to-M Reference Coordinates -------------------------------------
## Upstream reference minimum coverage values used to calculate 0% and 100%
## methylation reference ticks in M-space: M = log2((meth + pc) / (unmeth + pc))
min_CA_sites <- 3
sample_mvalue_reference_overrides <- list(
  "Chicken"        = list(CG_coverage = 5, CA_coverage = 5*(min_CA_sites), pseudocount = 0.1),
  "Frog"           = list(CG_coverage = 5, CA_coverage = 5*(min_CA_sites), pseudocount = 0.1),
  "Platypus"       = list(CG_coverage = 5, CA_coverage = 5*(min_CA_sites), pseudocount = 0.1),
  "Elephant Shark" = list(CG_coverage = 5, CA_coverage = 5*(min_CA_sites), pseudocount = 0.1),
  "Zebrafish"      = list(CG_coverage = 5, CA_coverage = 5*(min_CA_sites), pseudocount = 0.1)
)

## Reference overlay tick styling (short visible ticks near axes)
overlay_beta_reference_ticks               <- TRUE    # Draw short inset ticks at beta = 0% and 100%
beta_reference_overlay_tick_colour         <- "black" # Tick mark colour
beta_reference_overlay_tick_linewidth      <- 0.8     # Tick mark line thickness
beta_reference_overlay_tick_linetype       <- "solid" # Tick mark line style
beta_reference_overlay_tick_inset_fraction <- 0.001   # Distance of tick base from inner panel edge
beta_reference_overlay_tick_fraction       <- 0.01    # Tick mark length as fraction of axis span

## Reference label styling ("B-val = 0" and "B-val = 100")
show_beta_reference_labels                 <- TRUE    # Display textual reference labels
beta_reference_label_size                  <- 4       # Text font size
beta_reference_label_inset_fraction        <- 0.030   # Inset distance from axis edge into panel
beta_reference_label_side_offset_fraction  <- 0.06    # Offset along tick axis to prevent collision

## ---- 3. Grid Geometry & Axis Scaling ----------------------------------------
mvalue_bin_width_x           <- 0.2    # Bin width along x-axis (mCG M-value)
mvalue_bin_width_y           <- 0.2    # Bin width along y-axis (mCA M-value)
mvalue_bin_origin_x          <- 0.0    # Zero-centered grid alignment origin (X)
mvalue_bin_origin_y          <- 0.0    # Zero-centered grid alignment origin (Y)

mvalue_limit_quantiles       <- c(0, 1) # Quantile bounds for automatic axis range (full data extent)
mvalue_axis_padding_fraction <- 0.03    # Padding fraction added beyond extreme bins
force_square_panel           <- TRUE    # Enforce 1:1 square aspect ratio for data panel
coordinate_ratio             <- 1.0     # Coordinate aspect ratio

## ---- 4. Tile Count Transformation & Palette ---------------------------------
tile_fill_scale         <- "log2"   # Count transformation: linear, log2, log10, sqrt, asinh, relative
fill_palette            <- "magma"  # Continuous viridis palette family (e.g., "magma")
fill_direction          <- -1L      # -1 for reversed palette (lighter at low counts, dark at high)
fill_colour_gamma       <- 0.3      # Pacing exponent for high-density contrast (< 1 emphasizes dense regions)

apply_gradual_white_mix <- TRUE     # Mix gradual white into low counts for contrast
white_mix_low_values    <- 0.2      # Fraction of white added at lowest counts
white_mix_high_values   <- 0.1      # Fraction of white added at highest counts
white_mix_curve         <- 1.0      # Linear curve (1.0) for white-mix gradient

## ---- 5. Count-Contour Options -----------------------------------------------
add_count_contours_mvalues <- TRUE   # Overlay count-density contour rings
contour_break_mode         <- "explicit" # Explicit breaks vs automatic steps
contour_breaks_explicit    <- c(7, 8, 9, 10, 11, 12, 13, 14) # Transformed count levels for contours
contour_colour             <- "white" # Contour line colour
contour_linewidth          <- 0.8     # Contour line thickness
contour_linetype           <- "solid" # Contour line style
contour_alpha              <- 0.6     # Contour line transparency (0.0 to 1.0)

contour_smoothing          <- "gaussian" # Smoothing algorithm for contours (none, gaussian)
contour_smoothing_sigma    <- 1.0        # Standard deviation of 2D Gaussian kernel
contour_smoothing_radius   <- 3L         # Kernel truncation radius
contour_smoothing_edge     <- "reflect"  # Boundary condition for convolution
contour_include_zero_bins  <- TRUE       # Include empty bins in density surface

## ---- 6. Linear Model & Output Toggles ---------------------------------------
fit_linear_model         <- TRUE    # Fit OLS regression line (mCA ~ mCG)
lm_line_width            <- 1.0     # Regression line thickness
lm_line_colour           <- "grey1" # Regression line colour
lm_line_type             <- "solid" # Regression line style
lm_use_only_plotted_data <- FALSE   # TRUE = restrict fit to visible axis limits, FALSE = all data
show_lm_annotation       <- FALSE   # Annotate equation, R², and p-value on plot
lm_annotation_size       <- 4.0     # Annotation text font size

show_fill_legend         <- TRUE    # Display continuous colorbar legend
plot_width_in            <- 9.0     # Exported figure width (inches)
plot_height_in           <- 9.0     # Exported figure height (inches)
dpi                      <- 300     # PNG raster resolution
export_full_lm_summaries <- TRUE    # Export comprehensive text summaries of regression models

## =============================================================================
## EXECUTION PIPELINE & HELPER FUNCTIONS ----
## =============================================================================

normalise <- function(x) tolower(trimws(as.character(x)[1]))
num_token <- function(x) {
  z <- format(x, trim = TRUE, scientific = FALSE)
  gsub("-", "m", gsub("\\.", "p", z))
}
safe_stem <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

validate_settings <- function() {
  bool_flags <- list(
    force_square_panel = force_square_panel,
    show_fill_legend = show_fill_legend,
    fit_linear_model = fit_linear_model,
    apply_gradual_white_mix = apply_gradual_white_mix,
    add_count_contours_mvalues = add_count_contours_mvalues,
    contour_include_zero_bins = contour_include_zero_bins,
    overlay_beta_reference_ticks = overlay_beta_reference_ticks,
    show_beta_reference_labels = show_beta_reference_labels,
    lm_use_only_plotted_data = lm_use_only_plotted_data,
    show_lm_annotation = show_lm_annotation,
    export_full_lm_summaries = export_full_lm_summaries
  )
  for (nm in names(bool_flags)) {
    val <- bool_flags[[nm]]
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop(sprintf("Setting '%s' must be TRUE or FALSE. Got: %s", nm, class(val)[1]), call. = FALSE)
    }
  }
  if (!is.numeric(c(plot_width_in, plot_height_in, dpi, coordinate_ratio,
                    white_mix_low_values, white_mix_high_values,
                    contour_linewidth, contour_alpha, lm_line_width,
                    lm_annotation_size, beta_reference_overlay_tick_linewidth)) ||
      any(!is.finite(c(plot_width_in, plot_height_in, dpi, coordinate_ratio,
                       white_mix_low_values, white_mix_high_values,
                       contour_linewidth, contour_alpha, lm_line_width,
                       lm_annotation_size, beta_reference_overlay_tick_linewidth)))) {
    stop("Numeric settings must be finite numbers.", call. = FALSE)
  }
  if (contour_linewidth <= 0 || contour_alpha < 0 || contour_alpha > 1) {
    stop("Contour linewidth must be positive and alpha between 0 and 1.", call. = FALSE)
  }
  if (!normalise(tile_fill_scale) %in%
      c("linear", "log2", "log10", "sqrt", "asinh", "relative")) {
    stop("Unsupported tile_fill_scale.", call. = FALSE)
  }
  positive <- c(mvalue_bin_width_x, mvalue_bin_width_y, fill_colour_gamma, white_mix_curve)
  if (any(!is.finite(positive)) || any(positive <= 0)) {
    stop("Bin widths, gamma, and white_mix_curve must be positive.", call. = FALSE)
  }
  if (length(mvalue_limit_quantiles) != 2L ||
      mvalue_limit_quantiles[1] >= mvalue_limit_quantiles[2]) {
    stop("mvalue_limit_quantiles must have two increasing values.", call. = FALSE)
  }
}

read_windows <- function(path) {
  dt <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    cmd <- if (nzchar(Sys.which("pigz"))) paste("pigz -dc", shQuote(path))
    else paste("gzip -dc", shQuote(path))
    fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = FALSE)
  } else fread(path, sep = "\t", header = TRUE, showProgress = FALSE)

  aliases <- c(
    seqname = "chr", CG_mvalue = "CG_mvalues", CA_mvalue = "CA_mvalues"
  )
  for (from in names(aliases)) {
    to <- aliases[[from]]
    if (from %in% names(dt) && !(to %in% names(dt))) setnames(dt, from, to)
  }
  needed <- c("CG_mvalues", "CA_mvalues")
  missing <- setdiff(needed, names(dt))
  if (length(missing)) stop("Missing columns in ", basename(path), ": ",
                            paste(missing, collapse = ", "), call. = FALSE)
  dt[, (needed) := lapply(.SD, as.numeric), .SDcols = needed]
  dt
}

label_from_path <- function(path) {
  bn <- sub("\\.gz$", "", basename(path), ignore.case = TRUE)
  if (bn %in% names(sample_label_overrides)) return(sample_label_overrides[[bn]])
  lab <- sub("\\.tsv$", "", bn, ignore.case = TRUE)
  sub("^cowindows_?", "", lab, ignore.case = TRUE)
}

plot_option_stem <- function() {
  bin_token <- paste0("mxbw", num_token(mvalue_bin_width_x), "_mybw", num_token(mvalue_bin_width_y))
  paste("Mvalues", bin_token, normalise(tile_fill_scale),
        fill_palette, paste0("gamma", num_token(fill_colour_gamma)), sep = "_")
}

## ---- B-value reference helpers ----------------------------------------------
beta_extreme_mvalues <- function(coverage, pseudocount) {
  if (!is.finite(coverage) || coverage <= 0 ||
      !is.finite(pseudocount) || pseudocount <= 0) {
    return(c(beta0 = NA_real_, beta100 = NA_real_))
  }
  extreme <- log2((coverage + pseudocount) / pseudocount)
  c(beta0 = -extreme, beta100 = extreme)
}

sample_beta_reference <- function(sample) {
  cfg <- sample_mvalue_reference_overrides[[sample]]
  if (is.null(cfg)) return(list(
    CG_coverage = NA_real_, CA_coverage = NA_real_, pseudocount = NA_real_,
    x = c(beta0 = NA_real_, beta100 = NA_real_),
    y = c(beta0 = NA_real_, beta100 = NA_real_)))
  pc <- as.numeric(cfg$pseudocount)
  cg <- as.numeric(cfg$CG_coverage)
  ca <- as.numeric(cfg$CA_coverage)
  list(CG_coverage = cg, CA_coverage = ca, pseudocount = pc,
       x = beta_extreme_mvalues(cg, pc),
       y = beta_extreme_mvalues(ca, pc))
}

add_beta_reference_overlay_ticks <- function(p, reference, xlim, ylim) {
  if (!isTRUE(overlay_beta_reference_ticks)) return(p)

  xspan <- diff(xlim)
  yspan <- diff(ylim)

  ## mCG ticks, fully inset from the lower panel boundary
  for (nm in names(reference$x)) {
    xpos <- unname(reference$x[[nm]])
    if (!is.finite(xpos) || xpos < xlim[1] || xpos > xlim[2]) next

    tick_y_start <- ylim[1] + beta_reference_overlay_tick_inset_fraction * yspan
    tick_y_end   <- tick_y_start + beta_reference_overlay_tick_fraction * yspan

    p <- p + annotate(
      "segment",
      x = xpos, xend = xpos,
      y = tick_y_start, yend = tick_y_end,
      colour = beta_reference_overlay_tick_colour,
      linewidth = beta_reference_overlay_tick_linewidth,
      linetype = beta_reference_overlay_tick_linetype,
      lineend = "butt"
    )
  }

  ## mCA ticks, fully inset from the left panel boundary
  for (nm in names(reference$y)) {
    ypos <- unname(reference$y[[nm]])
    if (!is.finite(ypos) || ypos < ylim[1] || ypos > ylim[2]) next

    tick_x_start <- xlim[1] + beta_reference_overlay_tick_inset_fraction * xspan
    tick_x_end   <- tick_x_start + beta_reference_overlay_tick_fraction * xspan

    p <- p + annotate(
      "segment",
      x = tick_x_start, xend = tick_x_end,
      y = ypos, yend = ypos,
      colour = beta_reference_overlay_tick_colour,
      linewidth = beta_reference_overlay_tick_linewidth,
      linetype = beta_reference_overlay_tick_linetype,
      lineend = "butt"
    )
  }
  p
}

add_beta_reference_labels <- function(p, reference, xlim, ylim) {
  if (!isTRUE(show_beta_reference_labels)) return(p)
  xspan <- diff(xlim)
  yspan <- diff(ylim)

  for (nm in names(reference$x)) {
    xpos <- unname(reference$x[[nm]])
    if (!is.finite(xpos) || xpos < xlim[1] || xpos > xlim[2]) next
    is_zero <- identical(nm, "beta0")
    x_label <- xpos + if (is_zero) {
      -beta_reference_label_side_offset_fraction * xspan
    } else {
       beta_reference_label_side_offset_fraction * xspan
    }
    p <- p + annotate(
      "text", x = x_label,
      y = ylim[1] + beta_reference_label_inset_fraction * yspan,
      label = if (is_zero) "B-val = 0" else "B-val = 100",
      angle = 0, hjust = if (is_zero) 1 else 0, vjust = 0,
      size = beta_reference_label_size,
      colour = "grey50"
    )
  }

  for (nm in names(reference$y)) {
    ypos <- unname(reference$y[[nm]])
    if (!is.finite(ypos) || ypos < ylim[1] || ypos > ylim[2]) next
    is_zero <- identical(nm, "beta0")
    y_label <- ypos + if (is_zero) {
      -beta_reference_label_side_offset_fraction * yspan
    } else {
       beta_reference_label_side_offset_fraction * yspan
    }
    p <- p + annotate(
      "text", x = xlim[1] + beta_reference_label_inset_fraction * xspan,
      y = y_label,
      label = if (is_zero) "B-val = 0" else "B-val = 100",
      angle = 90, hjust = if (is_zero) 1 else 0, vjust = 0,
      size = beta_reference_label_size,
      colour = "grey50"
    )
  }
  p
}

## ---- colour palette ---------------------------------------------------------
apply_white_mix <- function(cols) {
  if (!apply_gradual_white_mix) return(cols)
  n <- length(cols)
  pos <- if (n == 1L) 1 else seq(0, 1, length.out = n)
  amounts <- white_mix_high_values +
    (white_mix_low_values - white_mix_high_values) *
    (1 - pos) ^ white_mix_curve
  rgb <- col2rgb(cols) / 255
  mixed <- sweep(rgb, 2, 1 - amounts, `*`)
  mixed <- sweep(mixed, 2, amounts, `+`)
  grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ])
}

fill_scale <- function() {
  title <- switch(normalise(tile_fill_scale),
    linear = "n windows", log2 = "n windows\nlog2(count+1)",
    log10 = "n windows\nlog10(count+1)", sqrt = "n windows\nsqrt(count)",
    asinh = "n windows\nasinh(count)", relative = "relative density")

  n <- 256L
  source <- viridisLite::viridis(n, option = fill_palette,
                                 direction = if (fill_direction < 0) -1 else 1)
  q <- (seq_len(n) - 0.5) / n
  idx <- round((q ^ (1 / fill_colour_gamma)) * (n - 1L)) + 1L
  cols <- apply_white_mix(source[pmax(1L, pmin(n, idx))])
  scale_fill_gradientn(name = title, colours = cols,
                       guide = if (show_fill_legend) guide_colorbar() else "none")
}

## ---- bins -------------------------------------------------------------------
aligned_breaks <- function(limits, width, origin = 0) {
  lower <- origin + floor((limits[1] - origin) / width) * width
  upper <- origin + ceiling((limits[2] - origin) / width) * width
  br <- seq(lower, upper, by = width)
  if (tail(br, 1) < upper - sqrt(.Machine$double.eps)) br <- c(br, upper)
  br
}

build_tiles <- function(dt, xcol, ycol, x_breaks, y_breaks) {
  d <- dt[is.finite(get(xcol)) & is.finite(get(ycol)),
          .(x = get(xcol), y = get(ycol))]
  d <- d[x >= min(x_breaks) & x <= max(x_breaks) &
         y >= min(y_breaks) & y <= max(y_breaks)]
  if (!nrow(d)) return(data.table())
  d[, ix := findInterval(x, x_breaks, rightmost.closed = TRUE, all.inside = TRUE)]
  d[, iy := findInterval(y, y_breaks, rightmost.closed = TRUE, all.inside = TRUE)]
  ans <- d[, .(count = .N), by = .(ix, iy)]
  ans[, `:=`(xmin = x_breaks[ix], xmax = x_breaks[ix + 1L],
             ymin = y_breaks[iy], ymax = y_breaks[iy + 1L])]
  ans[, fill_val := switch(normalise(tile_fill_scale),
    linear = count, log2 = log2(count + 1), log10 = log10(count + 1),
    sqrt = sqrt(count), asinh = asinh(count), relative = count / max(count))]
  ans[]
}

## ---- count-contours ---------------------------------------------------------
gaussian_kernel_1d <- function(sigma, radius = contour_smoothing_radius) {
  half_width <- max(1L, ceiling(radius * sigma))
  positions <- seq(-half_width, half_width)
  kernel <- exp(-(positions^2) / (2 * sigma^2))
  kernel / sum(kernel)
}

reflect_indices <- function(indices, n) {
  if (n == 1L) return(rep(1L, length(indices)))
  period <- 2L * n - 2L
  idx <- (indices - 1L) %% period
  ifelse(idx < n, idx + 1L, period - idx + 1L)
}

convolve_rows_gaussian <- function(mat, kernel, edge = contour_smoothing_edge) {
  nr <- nrow(mat); nc <- ncol(mat); half_width <- (length(kernel) - 1L) %/% 2L
  offsets <- seq(-half_width, half_width)
  out <- matrix(0, nr, nc)
  for (col_idx in seq_len(nc)) {
    cols_idx <- col_idx + offsets
    cols_idx <- switch(normalise(edge),
      reflect = reflect_indices(cols_idx, nc),
      clamp = pmin(pmax(cols_idx, 1L), nc),
      pmin(pmax(cols_idx, 1L), nc))
    out[, col_idx] <- rowSums(mat[, cols_idx, drop = FALSE] *
                              rep(kernel, each = nr))
  }
  out
}

convolve_cols_gaussian <- function(mat, kernel, edge = contour_smoothing_edge) {
  nr <- nrow(mat); nc <- ncol(mat); half_width <- (length(kernel) - 1L) %/% 2L
  offsets <- seq(-half_width, half_width)
  out <- matrix(0, nr, nc)
  for (row_idx in seq_len(nr)) {
    rows_idx <- row_idx + offsets
    rows_idx <- switch(normalise(edge),
      reflect = reflect_indices(rows_idx, nr),
      clamp = pmin(pmax(rows_idx, 1L), nr),
      pmin(pmax(rows_idx, 1L), nr))
    out[row_idx, ] <- colSums(mat[rows_idx, , drop = FALSE] *
                              rep(kernel, times = nc))
  }
  out
}

gaussian_smooth_2d <- function(mat, sigma) {
  k <- gaussian_kernel_1d(sigma)
  convolve_cols_gaussian(convolve_rows_gaussian(mat, k), k)
}

tile_count_matrix <- function(tiles, x_breaks, y_breaks, include_zeros = TRUE) {
  nx <- length(x_breaks) - 1L; ny <- length(y_breaks) - 1L
  mat <- matrix(0, nrow = ny, ncol = nx)
  if (!nrow(tiles)) return(mat)
  x_mids <- (head(x_breaks, -1L) + tail(x_breaks, -1L)) / 2
  y_mids <- (head(y_breaks, -1L) + tail(y_breaks, -1L)) / 2
  ix <- match(round(tiles$xmin + (tiles$xmax - tiles$xmin)/2, 6), round(x_mids, 6))
  iy <- match(round(tiles$ymin + (tiles$ymax - tiles$ymin)/2, 6), round(y_mids, 6))
  valid <- !is.na(ix) & !is.na(iy)
  mat[cbind(iy[valid], ix[valid])] <- tiles$count[valid]
  mat
}

add_count_contours <- function(p, tiles, x_breaks, y_breaks, enabled) {
  if (!isTRUE(enabled) || !nrow(tiles)) return(p)
  mat <- tile_count_matrix(tiles, x_breaks, y_breaks, include_zeros = contour_include_zero_bins)
  smooth_mat <- switch(normalise(contour_smoothing),
    gaussian = gaussian_smooth_2d(mat, contour_smoothing_sigma),
    mat)
  transformed <- switch(normalise(tile_fill_scale),
    linear = smooth_mat, log2 = log2(smooth_mat + 1), log10 = log10(smooth_mat + 1),
    sqrt = sqrt(smooth_mat), asinh = asinh(smooth_mat), relative = smooth_mat / max(smooth_mat))
  x_mids <- (head(x_breaks, -1L) + tail(x_breaks, -1L)) / 2
  y_mids <- (head(y_breaks, -1L) + tail(y_breaks, -1L)) / 2
  grid_df <- expand.grid(x = x_mids, y = y_mids)
  grid_df$z <- as.vector(t(transformed))

  breaks <- contour_breaks_explicit
  p + geom_contour(data = grid_df, aes(x = x, y = y, z = z),
                   breaks = breaks, colour = contour_colour,
                   linewidth = contour_linewidth, linetype = contour_linetype,
                   alpha = contour_alpha, inherit.aes = FALSE)
}

## ---- models -----------------------------------------------------------------
fit_model <- function(dt, xcol, ycol, sample, value_type,
                      xlim, ylim, bin_w_x, bin_w_y, n_bins_x, n_bins_y) {
  d <- dt[is.finite(get(xcol)) & is.finite(get(ycol)),
          .(x = get(xcol), y = get(ycol))]
  if (lm_use_only_plotted_data) {
    d <- d[x >= xlim[1] & x <= xlim[2] & y >= ylim[1] & y <= ylim[2]]
  }
  if (nrow(d) < 2L || uniqueN(d$x) < 2L) {
    return(list(model = NULL, stats = NULL))
  }
  fit <- lm(y ~ x, data = d)
  cf <- coef(fit); sm <- summary(fit)
  fstat <- sm$fstatistic
  pval <- if (!is.null(fstat)) {
    pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  } else NA_real_

  st <- data.table(
    sample = sample, value_type = value_type, n_points = nrow(d),
    intercept = unname(cf[1]), slope = unname(cf[2]),
    r_squared = sm$r.squared, adj_r_squared = sm$adj.r.squared,
    p_value = pval, sigma = sm$sigma,
    x_min = xlim[1], x_max = xlim[2], y_min = ylim[1], y_max = ylim[2],
    bin_width_x = bin_w_x, bin_width_y = bin_w_y,
    n_bins_x = n_bins_x, n_bins_y = n_bins_y
  )
  list(model = fit, stats = st)
}

add_model <- function(p, res, xlim, ylim) {
  if (is.null(res$model)) return(p)
  cf <- coef(res$model)
  p <- p + geom_abline(intercept = cf[1], slope = cf[2],
                       colour = lm_line_colour, linewidth = lm_line_width,
                       linetype = lm_line_type)
  if (show_lm_annotation && !is.null(res$stats)) {
    st <- res$stats
    txt <- sprintf("y = %.3f + %.3fx\nR² = %.3f; p = %.2e",
                   st$intercept, st$slope, st$r_squared, st$p_value)
    p <- p + annotate("text", x = xlim[1] + 0.05 * diff(xlim),
                      y = ylim[2] - 0.08 * diff(ylim),
                      label = txt, hjust = 0, vjust = 1,
                      size = lm_annotation_size, colour = lm_line_colour)
  }
  p
}

## ---- plotting ---------------------------------------------------------------
base_theme <- function() theme_classic(base_size = 22) + theme(
  plot.title = element_text(hjust = 0.5, size = 30),
  axis.text = element_text(size = 22, colour = "black"),
  axis.title = element_text(size = 28),
  legend.position = if (show_fill_legend) "right" else "none"
)

add_panel_coordinates <- function(p, xlim, ylim, force_square, fixed_ratio) {
  if (isTRUE(force_square)) {
    return(p + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
             theme(aspect.ratio = 1))
  }
  p + coord_fixed(ratio = fixed_ratio, xlim = xlim, ylim = ylim,
                  expand = FALSE, clip = "on")
}

build_m_plot <- function(dt, sample, xb, yb, model_result, beta_reference) {
  tiles <- build_tiles(dt, "CG_mvalues", "CA_mvalues", xb, yb)
  xlim <- range(xb)
  ylim <- range(yb)

  ## Layer order:
  ##   1. Opaque tile count density data
  ##   2. Count contours (optional)
  ##   3. Linear regression line (optional)
  ##   4. Inset overlay reference ticks (beta = 0% and 100%)
  ##   5. Inset reference labels ("B-val = 0" and "B-val = 100")
  p <- ggplot() +
    geom_rect(
      data = tiles,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_val),
      colour = NA, alpha = 1, inherit.aes = FALSE
    ) +
    fill_scale() +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0))

  p <- add_count_contours(
    p, tiles, xb, yb, enabled = add_count_contours_mvalues
  )
  if (fit_linear_model) p <- add_model(p, model_result, xlim, ylim)

  p <- add_beta_reference_overlay_ticks(p, beta_reference, xlim, ylim)
  p <- add_beta_reference_labels(p, beta_reference, xlim, ylim)

  p <- p +
    base_theme() +
    labs(title = paste0(sample), x = "mCG (M-value)", y = "mCA (M-value)")
  p <- add_panel_coordinates(
    p, xlim, ylim, force_square_panel, coordinate_ratio
  )
  p
}

## ---- Main Execution ---------------------------------------------------------
validate_settings()
paths <- list.files(input_dir, "^cowindows_.*\\.tsv(\\.gz)?$",
                    full.names = TRUE, ignore.case = TRUE)
if (!length(paths)) stop("No cowindows TSV files found in ", input_dir)

samples <- lapply(paths, function(path) list(
  path = path, sample = label_from_path(path), dt = read_windows(path)))

all_xm <- unlist(lapply(samples, function(s) s$dt$CG_mvalues), use.names = FALSE)
all_ym <- unlist(lapply(samples, function(s) s$dt$CA_mvalues), use.names = FALSE)
all_xm <- all_xm[is.finite(all_xm)]; all_ym <- all_ym[is.finite(all_ym)]
raw_xlim <- as.numeric(quantile(all_xm, mvalue_limit_quantiles, na.rm = TRUE))
raw_ylim <- as.numeric(quantile(all_ym, mvalue_limit_quantiles, na.rm = TRUE))
pad <- function(z) z + c(-1, 1) * max(diff(z), 1e-9) * mvalue_axis_padding_fraction
raw_xlim <- pad(raw_xlim); raw_ylim <- pad(raw_ylim)
m_xbreaks <- aligned_breaks(raw_xlim, mvalue_bin_width_x, mvalue_bin_origin_x)
m_ybreaks <- aligned_breaks(raw_ylim, mvalue_bin_width_y, mvalue_bin_origin_y)
m_xlim <- range(m_xbreaks); m_ylim <- range(m_ybreaks)

stats_list <- list(); full_models <- list(); reference_stats <- list()
for (i in seq_along(samples)) {
  s <- samples[[i]]; sample_stem <- safe_stem(s$sample)
  beta_reference <- sample_beta_reference(s$sample)
  reference_stats[[length(reference_stats) + 1L]] <- data.table(
    sample = s$sample,
    CG_reference_coverage = beta_reference$CG_coverage,
    CA_reference_coverage = beta_reference$CA_coverage,
    mvalue_pseudocount = beta_reference$pseudocount,
    CG_M_at_Bval_0 = unname(beta_reference$x["beta0"]),
    CG_M_at_Bval_100 = unname(beta_reference$x["beta100"]),
    CA_M_at_Bval_0 = unname(beta_reference$y["beta0"]),
    CA_M_at_Bval_100 = unname(beta_reference$y["beta100"])
  )

  mmodel <- fit_model(s$dt, "CG_mvalues", "CA_mvalues", s$sample, "M-values",
    m_xlim, m_ylim, mvalue_bin_width_x, mvalue_bin_width_y,
    length(m_xbreaks) - 1L, length(m_ybreaks) - 1L)

  if (!is.null(mmodel$stats)) stats_list[[length(stats_list) + 1L]] <- mmodel$stats
  full_models[[paste0(sample_stem, "_Mvalues")]] <- mmodel$model

  p <- build_m_plot(s$dt, s$sample, m_xbreaks, m_ybreaks, mmodel, beta_reference)
  path_png <- file.path(out_dir, paste0(sample_stem, "_", plot_option_stem(), ".png"))
  path_pdf <- file.path(out_dir, paste0(sample_stem, "_", plot_option_stem(), ".pdf"))
  ggsave(path_png, p, width = plot_width_in, height = plot_height_in, dpi = dpi)
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(path_pdf, p, width = plot_width_in, height = plot_height_in, device = pdf_dev)
}

reference_dt <- rbindlist(reference_stats, fill = TRUE)
fwrite(reference_dt, file.path(out_dir, "Bvalue_to_Mvalue_reference_positions.tsv"),
       sep = "\t", na = "NA")

stats <- rbindlist(stats_list, fill = TRUE)
stats[, plot_option_stem := plot_option_stem()]
fwrite(stats, file.path(out_dir, paste0(
  "linear_model_stats_Mbin", num_token(mvalue_bin_width_x), "x",
  num_token(mvalue_bin_width_y), ".tsv")), sep = "\t", na = "NA")

if (export_full_lm_summaries) {
  path <- file.path(out_dir, paste0(
    "linear_model_full_summaries_Mbin", num_token(mvalue_bin_width_x), "x",
    num_token(mvalue_bin_width_y), ".txt"))
  con <- file(path, "wt")
  for (nm in names(full_models)) {
    writeLines(c(strrep("=", 80), nm, strrep("=", 80)), con)
    model <- full_models[[nm]]
    if (is.null(model)) writeLines("Model unavailable.\n", con)
    else writeLines(c(capture.output(summary(model)),
                      paste("AIC:", AIC(model)), paste("BIC:", BIC(model)), ""), con)
  }
  close(con)
}

message("Done. M-value plots exported to: ", out_dir)
