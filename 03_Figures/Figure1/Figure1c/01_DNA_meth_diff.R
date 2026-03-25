# ================================
# Libraries
# ================================
library(readr)
library(GenomicRanges)

# ================================
# Helper functions
# ================================

# Read file and convert to GRanges
read_granges <- function(path) {
  df <- read_delim(path, delim = "\t")
  makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
}

# Merge overlapping CA/CG regions
merge_overlaps <- function(cg_df, ca_df, cg_gr, ca_gr) {
  overlaps <- findOverlaps(ca_gr, cg_gr)
  combined <- cbind(cg_df[subjectHits(overlaps), ], ca_df[queryHits(overlaps), ])
  combined$region_id <- paste(gsub("chr", "", combined$chr), combined$start, combined$end, sep = "-")
  row.names(combined) <- combined$region_id
  combined[, -c(11:13, 21)]
}

# Calculate mean methylation per sample for given columns
calculate_means <- function(df, cols) {
  colMeans(df[, cols], na.rm = TRUE) * 100  # multiply by 100
}

# Extract categories based on thresholds
get_categories <- function(df, cols, type = "cg") {
  if(type == "cg") {
    list(
      high = df[apply(df[, cols], 1, function(x) all(x > 0.66)), ],
      intermediate = df[apply(df[, cols], 1, function(x) all(x >= 0.33 & x <= 0.66)), ],
      low = df[apply(df[, cols], 1, function(x) all(x < 0.33)), ]
    )
  } else {
    df$mean <- rowMeans(df[, cols])
    q3 <- quantile(df$mean, 0.33)
    q1 <- quantile(df$mean, 0.66)
    list(
      high = df[apply(df[, cols], 1, function(x) all(x > q1)), ],
      intermediate = df[apply(df[, cols], 1, function(x) all(x >= q3 & x <= q1)), ],
      low = df[apply(df[, cols], 1, function(x) all(x < q3)), ]
    )
  }
}

# Difference between young and old
diff_young_old <- function(x) {
  young <- mean(x[1:3])
  old <- mean(x[5:7])
  old - young
}

# ================================
# CG context analysis (mCA depending on mCG)
# ================================
cg_bval <- "/home/epi-lab/Escritorio/Resultados/ventanas_cg_comun_cov5_bvalues"
ca_bval <- "/home/epi-lab/Escritorio/Resultados/ventanas_ca_comun_cov5_bvalues"

cg <- read_delim(cg_bval, delim = "\t")
ca <- read_delim(ca_bval, delim = "\t")
cg_gr <- makeGRangesFromDataFrame(cg, keep.extra.columns = TRUE)
ca_gr <- makeGRangesFromDataFrame(ca, keep.extra.columns = TRUE)

combined_windows <- merge_overlaps(cg, ca, cg_gr, ca_gr)
cg_categories <- get_categories(combined_windows, 4:10, type = "cg")

medias_mca_mcg <- calculate_means(cg_categories$high, 11:17)
medias_mca_hmcg <- calculate_means(cg_categories$intermediate, 11:17)
medias_mca_umcg <- calculate_means(cg_categories$low, 11:17)

diff_mca_mcg <- diff_young_old(medias_mca_mcg)
diff_mca_hmcg <- diff_young_old(medias_mca_hmcg)
diff_mca_umcg <- diff_young_old(medias_mca_umcg)

# ================================
# CA context analysis (mCG depending on mCA)
# ================================
cg <- read_delim(cg_bval, delim = "\t")
ca <- read_delim(ca_bval, delim = "\t")
cg_gr <- makeGRangesFromDataFrame(cg, keep.extra.columns = TRUE)
ca_gr <- makeGRangesFromDataFrame(ca, keep.extra.columns = TRUE)

combined_windows <- merge_overlaps(cg, ca, cg_gr, ca_gr)
ca_categories <- get_categories(combined_windows, 11:17, type = "ca")

medias_mcg_mca <- calculate_means(ca_categories$high, 4:10)
medias_mcg_hmca <- calculate_means(ca_categories$intermediate, 4:10)
medias_mcg_umca <- calculate_means(ca_categories$low, 4:10)

diff_mcg_mca <- diff_young_old(medias_mcg_mca)
diff_mcg_hmca <- diff_young_old(medias_mcg_hmca)
diff_mcg_umca <- diff_young_old(medias_mcg_umca)

# ================================
# Create summary dataframes
# ================================
mCA_context <- data.frame(
  bvalues = "mCA",
  umCG = diff_mca_umcg,
  hmCG = diff_mca_hmcg,
  mCG = diff_mca_mcg
)

mCG_context <- data.frame(
  bvalues = "mCG",
  umCA = diff_mcg_umca,
  hmCA = diff_mcg_hmca,
  mCA = diff_mcg_mca
)

# Print results
print(mCA_context)
print(mCG_context)
