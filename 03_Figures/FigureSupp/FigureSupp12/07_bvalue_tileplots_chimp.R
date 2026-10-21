#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## Methylation Analysis: B-value Tile Plots (chimp)
##
## Description:
##   Generates 2D density tile plots of genomic methylation
##   beta values (mCG vs mCA) along with a binned mean mCA trendline
##   scaled against a secondary right-hand y-axis.
##
## Environment:
##  Executed in R version 4.6.1.
##  Ubuntu 20.04.6 LTS
##
## Code Structure:
##   1. USER CONFIGURATION (Top Section):
##      All user-tunable parameters, including input/output directories, sample
##      label overrides, grid binning geometry, colour scale settings, mean
##      trendline appearance, and export dimensions.
##   2. EXECUTION PIPELINE & HELPER FUNCTIONS (Bottom Section):
##      Modular functions for data loading, density tile construction, binned
##      mean calculation, ggplot2 visualization, and dual-format (PNG + Cairo PDF)
##      figure export.
## -----------------------------------------------------------------------------

library(data.table)
library(ggplot2)

## =============================================================================
## USER CONFIGURATION ----
## =============================================================================

## ---- 1. File Paths & Dataset Settings ---------------------------------------
## Input directory containing per-sample mCG + mCA window files ('cowindows_*.tsv')
input_dir <- file.path("path/to/data/directory", "windows_chimpanzee", "resultados_por_muestra")

## Output directory for exported figures (PNG and PDF)
out_dir <- file.path("path/to/output/directory", "bvalue_tileplots_chimp")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Human-readable label overrides for sample plot titles
sample_label_overrides <- c(
  "cowindows_35y.tsv" = "Chimpanzee\n35 years old",
  "cowindows_40y.tsv" = "Chimpanzee\n40 years old",
  "cowindows_43y.tsv" = "Chimpanzee\n43 years old",
  "cowindows_44y.tsv" = "Chimpanzee\n44 years old",
  "cowindows_45y.tsv" = "Chimpanzee\n45 years old",
  "cowindows_57y.tsv" = "Chimpanzee\n57 years old"
)

## ---- 2. Grid Geometry & Axis Limits -----------------------------------------
bvalue_bin_width_x            <- 2        # Bin width along x-axis (mCG beta-value, percentage units)
bvalue_bin_width_y            <- 2        # Bin width along y-axis (mCA beta-value, percentage units)
bvalue_x_limits               <- c(0, 100) # Full axis limits for mCG (0% to 100%)
bvalue_y_limits               <- c(0, 100) # Full axis limits for mCA (0% to 100%)
force_square_panel            <- FALSE    # Enforce square coordinates via coord_cartesian aspect ratio
coordinate_ratio              <- 1        # Aspect ratio for fixed panel coordinates

## ---- 3. Tile Count Transformation & Palette ---------------------------------
fill_palette                  <- "magma"  # Continuous viridis palette family (e.g., "magma")
fill_direction                <- -1L      # -1 for reversed palette (lighter at low counts, dark at high)
fill_colour_gamma             <- 0.8      # Pacing exponent (< 1.0 enhances contrast in high-density regions)

apply_gradual_white_mix       <- TRUE     # Blend gradual white into low counts for background contrast
white_mix_low_values          <- 0.2      # Fraction of white blended into lowest count bins
white_mix_high_values         <- 0.1      # Fraction of white blended into highest count bins
white_mix_curve               <- 1.0      # Linear exponent for white mix gradient interpolation

## ---- 4. Binned Mean Trendline & Secondary Axis ------------------------------
show_mean_trendline           <- TRUE     # Compute and display binned mean trendline
mean_line_colour              <- "#0793e4" # Line and secondary axis colour (vibrant blue)
mean_line_width               <- 3.0      # Line thickness for the mean trendline
trend_bin_width               <- 2.0      # Bin width along mCG for calculating mean mCA values
remove_zero_secondary_y_label <- TRUE     # Suppress overlapping "0" label on secondary y-axis

## ---- 5. Output Figure Dimensions --------------------------------------------
show_fill_legend              <- TRUE     # Display the continuous count colourbar legend
plot_width_in                 <- 9.0      # Figure width in inches
plot_height_in                <- 9.0      # Figure height in inches
dpi                           <- 300      # Raster resolution (dots per inch) for PNG export

## =============================================================================
## EXECUTION PIPELINE & HELPER FUNCTIONS ----
## =============================================================================

num_token <- function(x) {
  z <- format(x, trim = TRUE, scientific = FALSE)
  gsub("-", "m", gsub("\\.", "p", z))
}
safe_stem <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

validate_settings <- function() {
  bool_flags <- list(
    force_square_panel = force_square_panel,
    show_fill_legend = show_fill_legend,
    apply_gradual_white_mix = apply_gradual_white_mix,
    show_mean_trendline = show_mean_trendline,
    remove_zero_secondary_y_label = remove_zero_secondary_y_label
  )
  for (nm in names(bool_flags)) {
    val <- bool_flags[[nm]]
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop(sprintf("Setting '%s' must be TRUE or FALSE. Got: %s", nm, class(val)[1]), call. = FALSE)
    }
  }
  if (!is.numeric(c(plot_width_in, plot_height_in, dpi, coordinate_ratio,
                    white_mix_low_values, white_mix_high_values,
                    mean_line_width, trend_bin_width)) ||
      any(!is.finite(c(plot_width_in, plot_height_in, dpi, coordinate_ratio,
                       white_mix_low_values, white_mix_high_values,
                       mean_line_width, trend_bin_width)))) {
    stop("Numeric settings must be finite numbers.", call. = FALSE)
  }
  positive <- c(bvalue_bin_width_x, bvalue_bin_width_y, fill_colour_gamma,
                white_mix_curve, mean_line_width, trend_bin_width)
  if (any(!is.finite(positive)) || any(positive <= 0)) {
    stop("Bin widths, gamma, mean_line_width, and trend_bin_width must be positive.", call. = FALSE)
  }
  if (length(bvalue_x_limits) != 2L || bvalue_x_limits[1] >= bvalue_x_limits[2] ||
      length(bvalue_y_limits) != 2L || bvalue_y_limits[1] >= bvalue_y_limits[2]) {
    stop("Axis limits must be length-2 numeric vectors with min < max.", call. = FALSE)
  }
}

read_windows <- function(path) {
  dt <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    cmd <- if (nzchar(Sys.which("pigz"))) paste("pigz -dc", shQuote(path))
    else paste("gzip -dc", shQuote(path))
    fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = FALSE)
  } else fread(path, sep = "\t", header = TRUE, showProgress = FALSE)

  aliases <- c(
    seqname = "chr", CG_bvalue = "CG_bvalues", CA_bvalue = "CA_bvalues"
  )
  for (from in names(aliases)) {
    to <- aliases[[from]]
    if (from %in% names(dt) && !(to %in% names(dt))) setnames(dt, from, to)
  }
  needed <- c("CG_bvalues", "CA_bvalues")
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
  bin_token <- paste0("bxbw", num_token(bvalue_bin_width_x), "_bybw", num_token(bvalue_bin_width_y))
  paste("Bvalues", bin_token, "log2", fill_palette,
        paste0("gamma", num_token(fill_colour_gamma)), sep = "_")
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
  title <- "n windows\nlog2(count+1)"
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
  ans[, fill_val := log2(count + 1)]
  ans[]
}

## ---- B-value mean trend -----------------------------------------------------
mean_trend <- function(dt) {
  d <- dt[is.finite(CG_bvalues) & is.finite(CA_bvalues)]
  d[, bin := cut(CG_bvalues, seq(0, 100, by = trend_bin_width), include.lowest = TRUE)]
  tr <- d[!is.na(bin), .(Mean = mean(CA_bvalues)), by = bin]
  tr[, x := as.numeric(sub("^[[(]([0-9.]+),.*$", "\\1", as.character(bin))) +
          trend_bin_width / 2]
  tr[order(x)]
}

## ---- theme and coordinates --------------------------------------------------
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

build_b_plot <- function(dt, sample, trend, trend_max) {
  xb <- aligned_breaks(bvalue_x_limits, bvalue_bin_width_x, 0)
  yb <- aligned_breaks(bvalue_y_limits, bvalue_bin_width_y, 0)
  tiles <- build_tiles(dt, "CG_bvalues", "CA_bvalues", xb, yb)

  sec_lab <- function(x) {
    z <- if (requireNamespace("scales", quietly = TRUE)) {
      scales::label_number()(x)
    } else {
      format(x, trim = TRUE)
    }
    if (remove_zero_secondary_y_label) z[abs(x) < 1e-10] <- ""
    z
  }

  p <- ggplot() +
    geom_rect(
      data = tiles,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_val),
      colour = NA, alpha = 1, inherit.aes = FALSE
    ) +
    fill_scale() +
    scale_x_continuous(
      limits = bvalue_x_limits, expand = c(0, 0),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_y_continuous(
      limits = bvalue_y_limits, expand = c(0, 0),
      breaks = c(0, 25, 50, 75, 100),
      sec.axis = sec_axis(~ . / 100 * trend_max,
                          name = "Mean mCA (β-value)",
                          labels = sec_lab)
    )

  if (show_mean_trendline && nrow(trend)) {
    tr <- copy(trend)[, yplot := Mean / trend_max * 100]
    p <- p + geom_line(data = tr, aes(x, yplot), colour = mean_line_colour,
                       linewidth = mean_line_width, inherit.aes = FALSE)
  }

  p <- p +
    base_theme() +
    theme(
      axis.text.y.right  = element_text(colour = mean_line_colour),
      axis.title.y.right = element_text(colour = mean_line_colour),
      axis.line.y.right  = element_line(colour = mean_line_colour),
      axis.ticks.y.right = element_line(colour = mean_line_colour)
    ) +
    labs(title = paste0(sample), x = "mCG (β-value)", y = "mCA (β-value)")

  p <- add_panel_coordinates(p, bvalue_x_limits, bvalue_y_limits,
                             force_square_panel, coordinate_ratio)
  p
}

## ---- Main Execution ---------------------------------------------------------
validate_settings()
paths <- list.files(input_dir, "^cowindows_.*\\.tsv(\\.gz)?$",
                    full.names = TRUE, ignore.case = TRUE)
if (!length(paths)) stop("No cowindows TSV files found in ", input_dir)

samples <- lapply(paths, function(path) list(
  path = path, sample = label_from_path(path), dt = read_windows(path)))

trends <- lapply(samples, function(s) mean_trend(s$dt))
trend_values <- unlist(lapply(trends, function(z) z$Mean), use.names = FALSE)
trend_max <- if (any(is.finite(trend_values))) max(trend_values, na.rm = TRUE) + 1 else 1

for (i in seq_along(samples)) {
  s <- samples[[i]]; sample_stem <- safe_stem(s$sample)
  p <- build_b_plot(s$dt, s$sample, trends[[i]], trend_max)
  path_png <- file.path(out_dir, paste0(sample_stem, "_", plot_option_stem(), ".png"))
  path_pdf <- file.path(out_dir, paste0(sample_stem, "_", plot_option_stem(), ".pdf"))
  ggsave(path_png, p, width = plot_width_in, height = plot_height_in, dpi = dpi)
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(path_pdf, p, width = plot_width_in, height = plot_height_in, device = pdf_dev)
}

message("Done. B-value plots exported to: ", out_dir)
