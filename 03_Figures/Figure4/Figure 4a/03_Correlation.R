library(readr)
library(dplyr)
library(GenomicRanges)

BASE_DIR <- "/home/epi-lab/Escritorio/Resultados/CHH/Correcciones/Size_windows/"
SAMPLES  <- c("35y", "43y", "55y", "67y", "78y", "85y", "103y")

OUT_DIR <- file.path(BASE_DIR, "resultados_por_muestra")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Common CA regions, same for all samples (loaded once)
comun_ca_GR <- makeGRangesFromDataFrame(
  read_delim(file.path(BASE_DIR, "CA_comun_1000"),
             delim = "\t", trim_ws = TRUE, show_col_types = FALSE),
  keep.extra.columns = TRUE
)

# Final summary table
summary_table <- data.frame()

for (sample in SAMPLES) {
  
  message("Processing: ", sample)
  
  # ── 1. Load CG ──────────────────────────────────────────────────────────────
  cpg_raw <- read_delim(
    file.path(BASE_DIR, paste0(  "CG_",sample,"_1000")),
    delim = "\t", trim_ws = TRUE, show_col_types = FALSE
  )
  cpg_GR <- makeGRangesFromDataFrame(cpg_raw, keep.extra.columns = TRUE)
  
  # Common CG regions for this sample
  comun_cpg_GR <- makeGRangesFromDataFrame(
    read_delim(
      file.path(BASE_DIR, paste0("CG_comun_1000bp.txt")),  # adjust name if needed
      delim = "\t", trim_ws = TRUE, show_col_types = FALSE
    ),
    keep.extra.columns = TRUE
  )
  
  # Filter CG by its own common regions
  cgsample       <- as.data.frame(cpg_GR[cpg_GR %over% comun_cpg_GR])
  cgsample$width <- NULL
  colnames(cgsample)[1] <- "chr"
  
  # CG metrics
  cgsample$CG_bvalues <- (cgsample$numCs / cgsample$coverage) * 100
  cgsample$CG_mvalues     <- log2((cgsample$numCs + 0.1) / (cgsample$numTs + 0.1))
  
  # Save CG windows (before overlapping with CA)
  write_tsv(
    cgsample[, c("chr", "start", "end", "CG_bvalues", "CG_mvalues")],
    file.path(OUT_DIR, paste0("CG_windows_", sample, ".tsv"))
  )
  
  # ── 2. Load CA ──────────────────────────────────────────────────────────────
  ca_raw <- read_delim(
    file.path(BASE_DIR, paste0( "CA_",sample,"_1000")),  # adjust name if needed
    delim = "\t", trim_ws = TRUE, show_col_types = FALSE
  )
  ca_raw$CA_bvalues <- (ca_raw$numCs / ca_raw$coverage) * 100
  ca_raw$CA_mvalues     <- log2((ca_raw$numCs + 0.1) / (ca_raw$numTs + 0.1))
  
  ca_GR <- makeGRangesFromDataFrame(ca_raw, keep.extra.columns = TRUE)
  
  # Filter CA by the shared common regions
  ca_filt       <- as.data.frame(ca_GR[ca_GR %over% comun_ca_GR])
  ca_filt$width <- NULL
  colnames(ca_filt)[1] <- "chr"
  
  # Save CA windows (before overlapping with CG categories)
  write_tsv(
    ca_filt[, c("chr", "start", "end", "CA_bvalues", "CA_mvalues")],
    file.path(OUT_DIR, paste0("CA_windows_", sample, ".tsv"))
  )
  
  # ── 3. Classify CG windows ──────────────────────────────────────────────────
  unmethylated_GR   <- makeGRangesFromDataFrame(
    filter(cgsample, CG_bvalues < 33), keep.extra.columns = TRUE)
  hemimethylated_GR <- makeGRangesFromDataFrame(
    filter(cgsample, CG_bvalues >= 33 & CG_bvalues <= 66), keep.extra.columns = TRUE)
  methylated_GR     <- makeGRangesFromDataFrame(
    filter(cgsample, CG_bvalues > 66), keep.extra.columns = TRUE)
  
  ca_GR_filt <- makeGRangesFromDataFrame(ca_filt, keep.extra.columns = TRUE)
  
  # ── 4. Overlaps CA ↔ CG classes ─────────────────────────────────────────────
  make_windows <- function(ca_gr, cg_gr) {
    hits <- findOverlaps(ca_gr, cg_gr)
    df   <- cbind(
      as.data.frame(cg_gr[subjectHits(hits)])[, c("seqnames", "start", "end", "CG_bvalues", "CG_mvalues")],
      as.data.frame(ca_gr[queryHits(hits)])[,  c("CA_bvalues", "CA_mvalues")]
    )
    colnames(df)[1] <- "chr"
    df
  }
  
  UMwindows <- make_windows(ca_GR_filt, unmethylated_GR)
  HMwindows <- make_windows(ca_GR_filt, hemimethylated_GR)
  Mwindows  <- make_windows(ca_GR_filt, methylated_GR)
  
  # cowindows with sample name
  assign(paste0("cowindows_", sample), rbind(UMwindows, HMwindows, Mwindows))
  cowindows <- get(paste0("cowindows_", sample))
  
  # Save cowindows for this sample
  write_tsv(
    cowindows,
    file.path(OUT_DIR, paste0("cowindows_", sample, ".tsv"))
  )
  
  # ── 5. Statistics ────────────────────────────────────────────────────────────
  foldchange  <- mean(Mwindows$CA_bvalues) / mean(UMwindows$CA_bvalues)
  spear       <- cor.test(cowindows$CA_bvalues, cowindows$CG_bvalues, method = "spearman")
  reg         <- summary(lm(CA_mvalues ~ CG_mvalues, data = cowindows))
  
  summary_table <- rbind(summary_table, data.frame(
    Sample       = sample,
    CG_lvl       = mean(cgsample$CG_bvalues),
    CA_lvl       = mean(ca_filt$CA_bvalues),
    mCA_mCG      = mean(Mwindows$CA_bvalues),
    mCA_hmCG     = mean(HMwindows$CA_bvalues),
    mCA_umCG     = mean(UMwindows$CA_bvalues),
    delta_mCA    = mean(Mwindows$CA_bvalues) - mean(ca_filt$CA_bvalues),
    FC           = foldchange,
    Spearman     = spear$estimate,
    pvalue_spear = spear$p.value,
    Slope        = reg$coefficients[2, 1],
    R_adj        = reg$adj.r.squared,
    Intercept    = reg$coefficients[1, 1],
    pvalue_lm    = reg$coefficients[2, 4]
  ))
}

# Show and save the summary table
print(summary_table)
write_tsv(summary_table, file.path(BASE_DIR, "summary_CA_CG_1000.tsv"))