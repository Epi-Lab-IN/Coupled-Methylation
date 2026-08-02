###################################################################
# Find common regions across cg_M1..M7 and ca_M1..M7 files and build
# 4 wide-format output tables (chr, start, end, M1..M7):
#   1) CG_bvalues_common.tsv
#   2) CG_mvalues_common.tsv
#   3) CA_bvalues_common.tsv
#   4) CA_mvalues_common.tsv
#
# cg and ca intersections are computed SEPARATELY (cg vs cg, ca vs ca).
###################################################################

## ----------------------------------------------------------------
## 0. CONFIG -- edit only this section
## ----------------------------------------------------------------

BASE_DIR <- "/ruta/a/tus/archivos/"   # <-- CHANGE this to your folder
SAMPLES  <- paste0("M", 1:7)          # M1, M2, ..., M7

OUT_DIR <- file.path(BASE_DIR, "resultados")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

library(data.table)

## ----------------------------------------------------------------
## 1. HELPER FUNCTIONS
## ----------------------------------------------------------------

# Reads cg_MX.txt / ca_MX.txt and normalizes "chr" to character
# (in case one file has it as text and another as numeric)
read_track <- function(prefix, sample, dir = BASE_DIR) {
  f <- file.path(dir, paste0(prefix, "_", sample, ".txt"))
  dt <- fread(f, header = TRUE, sep = "\t")
  dt[, chr := as.character(chr)]
  dt
}

# Builds a wide table (one column per sample) for a given value column,
# restricted to a set of common regions
build_wide <- function(dt_list, common_regions, value_col) {
  wide <- copy(common_regions)
  for (s in names(dt_list)) {
    matched <- dt_list[[s]][common_regions, on = .(chr, start, end)]
    wide[[s]] <- matched[[value_col]]
  }
  wide
}

## ----------------------------------------------------------------
## 2. Load all cg and ca files
## ----------------------------------------------------------------

cg_list <- lapply(SAMPLES, read_track, prefix = "cg")
ca_list <- lapply(SAMPLES, read_track, prefix = "ca")
names(cg_list) <- SAMPLES
names(ca_list) <- SAMPLES

## ----------------------------------------------------------------
## 3. Common regions within cg samples, and within ca samples
##    (cg and ca are NOT intersected with each other)
## ----------------------------------------------------------------

common_cg_regions <- Reduce(fintersect, lapply(cg_list, function(dt) unique(dt[, .(chr, start, end)])))
common_ca_regions <- Reduce(fintersect, lapply(ca_list, function(dt) unique(dt[, .(chr, start, end)])))

setorder(common_cg_regions, chr, start, end)
setorder(common_ca_regions, chr, start, end)

message("Common regions across cg samples: ", nrow(common_cg_regions))
message("Common regions across ca samples: ", nrow(common_ca_regions))

## ----------------------------------------------------------------
## 4. Build the 4 final wide-format tables
## ----------------------------------------------------------------

CG_bvalues <- build_wide(cg_list, common_cg_regions, "CG_bvalues")
CG_mvalues <- build_wide(cg_list, common_cg_regions, "CG_mvalues")
CA_bvalues <- build_wide(ca_list, common_ca_regions, "CA_bvalues")
CA_mvalues <- build_wide(ca_list, common_ca_regions, "CA_mvalues")

## ----------------------------------------------------------------
## 5. Save results
## ----------------------------------------------------------------

fwrite(CG_bvalues, file.path(OUT_DIR, "CG_bvalues_common.tsv"), sep = "\t")
fwrite(CG_mvalues, file.path(OUT_DIR, "CG_mvalues_common.tsv"), sep = "\t")
fwrite(CA_bvalues, file.path(OUT_DIR, "CA_bvalues_common.tsv"), sep = "\t")
fwrite(CA_mvalues, file.path(OUT_DIR, "CA_mvalues_common.tsv"), sep = "\t")

message("Done. Files saved in: ", OUT_DIR)