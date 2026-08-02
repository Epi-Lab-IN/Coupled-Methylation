# ==============================================================================
# Script: Cross-species liftOver and ChIP-seq enrichment analysis
# Description: Converts human (hg38) CoM+ regions to mouse (mm10) coordinates
#              using liftOver, then tests for enrichment of overlap between
#              the lifted regions and a ChIP-seq peak set using Fisher's
#              exact test.
# ==============================================================================

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# 1. Load human CoM regions (hg38) and convert to GRanges
# ------------------------------------------------------------------------------

com_pos_hg38 <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/COM_pos_filtrated",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

com_pos_hg38$chr <- paste0("chr", com_pos_hg38$chr)
com_pos_hg38_gr  <- makeGRangesFromDataFrame(com_pos_hg38, keep.extra.columns = TRUE)

# ------------------------------------------------------------------------------
# 2. LiftOver: convert CoM regions from human (hg38) to mouse (mm10)
# ------------------------------------------------------------------------------

chain_hg38_to_mm10 <- import.chain(
  "/home/epi-lab/Escritorio/Database/hg38ToMm10.over.chain"
)

com_pos_mm10_list <- liftOver(com_pos_hg38_gr, chain_hg38_to_mm10)
com_pos_mm10_gr    <- unlist(com_pos_mm10_list)  

write.table(
  as.data.frame(com_pos_mm10_gr),
  file = "CoM_pos_mouse_CHH",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------------------------
# 3. Helper function to load a BED-like file as a GRanges object
# ------------------------------------------------------------------------------

read_granges <- function(path) {
  df <- read_delim(path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
}

# ------------------------------------------------------------------------------
# 4. Build the background (universe) of mouse CpG/CpH positions
# ------------------------------------------------------------------------------

mouse_ca_common <- read_delim(
  "/home/epi-lab/Escritorio/Animals/CHH/Mouse/CAH_VC",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

mouse_cg_common <- read_delim(
  "/home/epi-lab/Escritorio/Animals/CHH/Mouse/CG_VC",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

mouse_background <- bind_rows(mouse_cg_common, mouse_ca_common) %>%
  distinct(chr, start, end, .keep_all = TRUE)

background_gr <- makeGRangesFromDataFrame(mouse_background, keep.extra.columns = TRUE)

# ------------------------------------------------------------------------------
# 5. Load lifted CoM target regions (mouse) and ChIP-seq peaks
# ------------------------------------------------------------------------------

target_mm10 <- read_delim(
  "/home/epi-lab/Escritorio/Resultados/CHH/Corr_age/CoM_pos_mouse_CHH",
  delim = "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE
)
target_gr <- makeGRangesFromDataFrame(target_mm10, keep.extra.columns = TRUE)

chip_peaks_gr <- read_granges(
  "/home/epi-lab/Escritorio/chip-seq/DNMT3A.txt"
)

# ------------------------------------------------------------------------------
# 6. Build 2x2 contingency table: overlap with ChIP peaks vs. target regions,
#    relative to the background universe
# ------------------------------------------------------------------------------

in_chip_in_target   <- sum((background_gr %over% chip_peaks_gr) &  (background_gr %over% target_gr))
in_chip_out_target  <- sum((background_gr %over% chip_peaks_gr) & !(background_gr %over% target_gr))
out_chip_in_target  <- sum(!(background_gr %over% chip_peaks_gr) &  (background_gr %over% target_gr))
out_chip_out_target <- sum(!(background_gr %over% chip_peaks_gr) & !(background_gr %over% target_gr))

contingency_table <- matrix(
  c(in_chip_in_target, in_chip_out_target,
    out_chip_in_target, out_chip_out_target),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    ChIP   = c("In ChIP", "Outside ChIP"),
    Target = c("In target", "Outside target")
  )
)

print(contingency_table)

# ------------------------------------------------------------------------------
# 7. Fisher's exact test: enrichment of ChIP peak overlap in target regions
# ------------------------------------------------------------------------------

fisher_result <- fisher.test(contingency_table)

cat("\nOdds ratio:", fisher_result$estimate, "\n")
cat("p-value:", fisher_result$p.value, "\n")