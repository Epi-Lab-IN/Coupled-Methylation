###################################################################
# Methylation ~ age correlation, per genomic window (CA and CG)
#
# Simple linear regression (value ~ age) computed with closed-form
# matrix formulas for ALL rows at once (no per-row lm() calls).
#
###################################################################

library(data.table)

## ----------------------------------------------------------------
## 0. CONFIG
## ----------------------------------------------------------------

ages_asc <- c(35, 43, 55, 67, 78, 85, 103)   # sample ages, in the column order of the file
P_CUTOFF <- 0.05

## ----------------------------------------------------------------
## 1. FUNCTION
## ----------------------------------------------------------------

run_age_correlation <- function(in_file, out_file, context, slope_direction,
                                ages = ages_asc, p_cutoff = P_CUTOFF) {
  
  suffix    <- tolower(context)                       # "ca" or "cg"
  col_names <- paste0("age", ages, "_", suffix)       # names in file order (columns 4:10)
  
  # Load data (fread is much faster than read_delim for large files)
  windows <- fread(in_file, sep = "\t")
  setnames(windows, 4:(3 + length(ages)), col_names)
  
  # Value matrix (rows = windows, cols = samples), z-scored per sample column
  Y <- scale(as.matrix(windows[, ..col_names]))
  
  # Vectorized simple linear regression for every row
  x   <- ages
  n   <- length(x)
  xc  <- x - mean(x)
  Sxx <- sum(xc^2)
  
  ybar <- rowMeans(Y)
  Sxy  <- as.vector(Y %*% xc)
  Syy  <- rowSums(Y^2) - n * ybar^2
  rm(Y)
  
  slope     <- Sxy / Sxx
  intercept <- ybar - slope * mean(x)
  SSE       <- pmax(Syy - slope * Sxy, 0)           # clip tiny negatives from floating point error
  se_slope  <- sqrt((SSE / (n - 2)) / Sxx)
  p_value   <- 2 * pt(-abs(slope / se_slope), df = n - 2)
  r_value   <- sign(slope) * sqrt(1 - SSE / Syy)
  
  results <- data.table(
    chr       = windows[[1]],
    start     = windows[[2]],
    end       = windows[[3]],
    R         = r_value,
    p.value   = p_value,
    slope     = slope,
    intercept = intercept
  )
  rm(windows)
  
  # Keep significant windows with the requested slope direction
  keep <- if (slope_direction == "positive") {
    results[p.value < p_cutoff & slope > 0]
  } else {
    results[p.value < p_cutoff & slope < 0]
  }
  
  fwrite(keep, out_file, sep = "\t")
  message(context, " (", slope_direction, " slope): ",
          nrow(keep), " significant windows out of ", nrow(results))
  
  invisible(keep)
}

## ----------------------------------------------------------------
## 2. RUN
## ----------------------------------------------------------------

res_ca <- run_age_correlation(
  in_file         = "/home/epi-lab/Escritorio/Resultados/CHH/CG_mvalues_common.tsv",  # check this path (see note)
  out_file        = "ca_correlation_positive.tsv",
  context         = "CA",
  slope_direction = "positive"
)

res_cg <- run_age_correlation(
  in_file         = "/home/epi-lab/Escritorio/Resultados/CHH/CG_mvalues_common.tsv",
  out_file        = "cg_correlation_negative.tsv",
  context         = "CG",
  slope_direction = "negative"
)