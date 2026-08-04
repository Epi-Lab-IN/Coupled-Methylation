## 02 — Portable window-cache and CoM-metric generation
##
## purpose: generate the BAM-downsampled metric table consumed by 03
## 03 is the plotting script.
##
## Required input layout (written by 01_downsample_bams.sh):
##   /path/to/01_results/<sample>/CX_processed/Data/*_CG.gz
##   /path/to/01_results/<sample>/CX_processed/Data/*_CA.gz
##
## Primary output:
##   /path/to/02_results/tables/all_samples_CoM_metrics_by_bam_depth_and_mincov.tsv

suppressPackageStartupMessages({
  library(data.table)
})

## ---- User configuration: edit paths before running --------------------------
context_root <- "/path/to/01_results"
outdir <- "/path/to/02_results"

# sample IDs and ages
sample_info <- data.table(
  sample = paste0("M", 1:7),
  age = c(35L, 43L, 55L, 67L, 78L, 85L, 103L),
  data_dir = file.path(context_root, paste0("M", 1:7), "CX_processed", "Data")
)

# configuration for window tiling and filtering

win_size <- 100L # window tiling width in base pairs
cg_cov_bases <- 1L # number of covered CG sites required to retain a window
ca_cov_bases <- 3L # number of covered CA sites required to retain a window
window_mincov_values <- c(1L, 5L, 10L) # minimum read coverage required at an individual cytosine site
autosomal_nums <- as.character(1:22) # chromosomes to include in the tiling and metrics
min_free_mem_gb <- 8 # minimum free RAM (GB) required to process a chromosome split file

## Set TRUE only when inputs or the settings above change. Otherwise valid cache
## files and the metric table are reused.
force_rebuild_chr_split <- TRUE
force_rebuild_windows <- TRUE
force_recompute_metrics <- TRUE

## ---- Output directories and run record --------------------------------------
table_dir <- file.path(outdir, "tables")
cache_dir <- file.path(outdir, "window_cache_directory")
split_root <- file.path(outdir, "chr_split_contextsafe")
log_dir <- file.path(outdir, "logs")
for (d in c(table_dir, cache_dir, split_root, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "02_generate_window_cache_and_metrics.log")
log_msg <- function(msg) {
  m <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", msg, "\n")
  cat(m); cat(m, file = log_file, append = TRUE)
}
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_start.txt"))
fwrite(data.table(
  parameter = c("context_root", "outdir", "win_size", "cg_cov_bases", "ca_cov_bases", "window_mincov_values", "force_rebuild_chr_split", "force_rebuild_windows", "force_recompute_metrics"),
  value = c(context_root, outdir, win_size, cg_cov_bases, ca_cov_bases, paste(window_mincov_values, collapse = ","), force_rebuild_chr_split, force_rebuild_windows, force_recompute_metrics)
), file.path(log_dir, "run_configuration.tsv"), sep = "\t")

## ---- Helpers ----------------------------------------------------------------
available_mem_gb <- function() {
  if (!file.exists("/proc/meminfo")) return(NA_real_)
  x <- grep("^MemAvailable:", readLines("/proc/meminfo", warn = FALSE), value = TRUE)
  if (length(x) == 1L) as.numeric(gsub("[^0-9]", "", x)) / 1024^2 else NA_real_
}
check_mem <- function(where = "") {
  x <- available_mem_gb()
  if (!is.na(x) && x < min_free_mem_gb) stop("Available RAM < ", min_free_mem_gb, " GB at ", where)
}
parse_depth_x <- function(path) suppressWarnings(as.numeric(sub("^.*_([0-9]+)x_.*$", "\\1", basename(path))))

split_context_file_by_chr <- function(input_file, context_label, sample_id, depth_x) {
  prefix <- paste0(sample_id, "_", depth_x, "x_", context_label)
  out_dir <- file.path(split_root, prefix)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  expected <- file.path(out_dir, paste0(prefix, "_chr_", autosomal_nums, ".txt.gz"))
  if (!isTRUE(force_rebuild_chr_split) && all(file.exists(expected)) && all(file.info(expected)$size > 0)) return(out_dir)

  log_msg(paste("Creating chromosome splits:", sample_id, depth_x, "x", context_label))
  unlink(file.path(out_dir, "*.txt")); unlink(file.path(out_dir, "*.txt.gz"))
  regex <- if (context_label == "CA") "^CA" else "^CG"
  awk_program <- 'BEGIN{OFS="\t"} {chr=$1; gsub(/^chr/,"",chr); if(chr>=1 && chr<=22 && $7 ~ re){outfile=dir "/" pref "_chr_" chr ".txt"; print >> outfile}}'
  cmd <- sprintf("zcat %s | awk -v dir=%s -v pref=%s -v re=%s %s",
                 shQuote(input_file), shQuote(out_dir), shQuote(prefix), shQuote(regex), shQuote(awk_program))
  if (!identical(system(cmd), 0L)) stop("Chromosome split failed: ", input_file)
  txt_files <- list.files(out_dir, pattern = "\\.txt$", full.names = TRUE)
  if (!length(txt_files)) stop("No chromosome-split files produced: ", input_file)
  if (!identical(system(paste("gzip -f", paste(shQuote(txt_files), collapse = " "))), 0L)) stop("gzip failed in ", out_dir)
  out_dir
}

tile_context_chr_split <- function(split_dir, context_label, sample_id, depth_x, mincov, cov_bases) {
  out <- file.path(cache_dir, paste0(sample_id, "_", depth_x, "x_", context_label,
                                     "_windows_", win_size, "bp_mincov", mincov, "_sites", cov_bases, ".tsv.gz"))
  if (file.exists(out) && !isTRUE(force_rebuild_windows) && file.info(out)$size > 0) return(out)
  prefix <- basename(split_dir)
  chr_files <- file.path(split_dir, paste0(prefix, "_chr_", autosomal_nums, ".txt.gz"))
  chr_files <- chr_files[file.exists(chr_files)]
  if (!length(chr_files)) stop("No chromosome-split files in ", split_dir)
  log_msg(paste("Tiling", sample_id, depth_x, "x", context_label, "mincov", mincov))

  pieces <- list()
  for (chr_file in chr_files) {
    check_mem(basename(chr_file))
    dt <- fread(chr_file, sep = "\t", header = FALSE, select = c(1, 2, 4, 5),
                col.names = c("chr", "pos", "count_M", "count_U"), showProgress = FALSE)
    if (!nrow(dt)) next
    dt[, `:=`(chr = as.character(chr), pos = as.integer(pos), count_M = as.numeric(count_M), count_U = as.numeric(count_U))]
    dt[, coverage_site := count_M + count_U]
    dt <- dt[coverage_site >= mincov]
    if (!nrow(dt)) { rm(dt); gc(); next }
    dt[, start := ((pos - 1L) %/% win_size) * win_size + 1L]
    dt[, end := start + win_size - 1L]
    agg <- dt[, .(coverage = sum(coverage_site), numCs = sum(count_M), numTs = sum(count_U), n_sites = .N), by = .(chr, start, end)]
    agg <- agg[n_sites >= cov_bases]
    if (nrow(agg)) pieces[[length(pieces) + 1L]] <- agg[, .(chr, start, end, strand = "*", coverage, numCs, numTs)]
    rm(dt, agg); gc()
  }
  if (length(pieces)) {
    ans <- rbindlist(pieces); setorder(ans, chr, start, end); fwrite(ans, out, sep = "\t")
  } else {
    fwrite(data.table(chr=character(), start=integer(), end=integer(), strand=character(), coverage=numeric(), numCs=numeric(), numTs=numeric()), out, sep="\t")
  }
  out
}

read_window_dt <- function(path) {
  x <- fread(path, showProgress = FALSE)
  if (nrow(x)) { x[, chr := as.character(chr)]; x[, `:=`(start=as.integer(start), end=as.integer(end))] }
  x
}

compute_com_metrics <- function(cg_file, ca_file, sample_id, age, depth_x, mincov) {
  cg <- read_window_dt(cg_file); ca <- read_window_dt(ca_file)
  base <- data.table(sample=sample_id, age=age, depth_x=depth_x, achieved_depth_x=depth_x,
                     mincov_cg=mincov, mincov_ca=mincov, pair_type="matched_min_coverage")
  if (!nrow(cg) || !nrow(ca)) return(cbind(base, data.table(n_windows=0L, slope=NA_real_, mCA_mean=NA_real_, mCG_mean=NA_real_)))
  cg[, `:=`(CG_methylation=100*numCs/coverage, mvalues_cg=log2((numCs+0.1)/(numTs+0.1)))]
  ca[, `:=`(CA_methylation=100*numCs/coverage, mvalues_ca=log2((numCs+0.1)/(numTs+0.1)))]
  common <- merge(cg[, .(chr,start,end,CG_methylation,mvalues_cg)], ca[, .(chr,start,end,CA_methylation,mvalues_ca)], by=c("chr","start","end"), all=FALSE, sort=FALSE)
  if (nrow(common) < 10L) return(cbind(base, data.table(n_windows=nrow(common), slope=NA_real_, mCA_mean=mean(common$CA_methylation, na.rm=TRUE), mCG_mean=mean(common$CG_methylation, na.rm=TRUE))))
  fit <- tryCatch(lm(mvalues_ca ~ mvalues_cg, data=common), error=function(e) NULL)
  cbind(base, data.table(n_windows=nrow(common), slope=if (is.null(fit)) NA_real_ else coef(fit)[2],
                          mCA_mean=mean(common$CA_methylation, na.rm=TRUE), mCG_mean=mean(common$CG_methylation, na.rm=TRUE)))
}

## ---- Discover paired CG/CA files --------------------------------------------
index_list <- list()
for (i in seq_len(nrow(sample_info))) {
  s <- sample_info$sample[i]; d <- sample_info$data_dir[i]
  if (!dir.exists(d)) stop("Missing context directory for ", s, ": ", d)
  cg <- list.files(d, pattern="_CG\\.gz$", full.names=TRUE)
  ca <- list.files(d, pattern="_CA\\.gz$", full.names=TRUE)
  cgdt <- data.table(sample=s, cg_file=cg, depth_x=parse_depth_x(cg))
  cadt <- data.table(sample=s, ca_file=ca, depth_x=parse_depth_x(ca))
  index_list[[length(index_list)+1L]] <- merge(cgdt, cadt, by=c("sample","depth_x"), all=FALSE)
}
file_index <- rbindlist(index_list, use.names=TRUE, fill=TRUE)
file_index <- merge(file_index, sample_info[, .(sample,age)], by="sample", all.x=TRUE)
file_index <- file_index[is.finite(depth_x)][order(sample, depth_x)]
if (!nrow(file_index)) stop("No paired CG/CA context files found.")
fwrite(file_index, file.path(table_dir, "all_samples_downsampled_bam_file_index.tsv"), sep="\t")

## ---- Split, tile, and calculate required metrics ----------------------------
metric_cache <- file.path(table_dir, "all_samples_CoM_metrics_by_bam_depth_and_mincov.tsv")
if (file.exists(metric_cache) && !isTRUE(force_recompute_metrics)) {
  metrics_dt <- fread(metric_cache)
  log_msg("Reusing existing metric table: ")
  log_msg(metric_cache)
} else {
  metrics_list <- list()
  for (i in seq_len(nrow(file_index))) {
    s <- file_index$sample[i]; depth <- file_index$depth_x[i]
    cg_split <- split_context_file_by_chr(file_index$cg_file[i], "CG", s, depth)
    ca_split <- split_context_file_by_chr(file_index$ca_file[i], "CA", s, depth)
    for (mincov in window_mincov_values) {
      cg_window <- tile_context_chr_split(cg_split, "CG", s, depth, mincov, cg_cov_bases)
      ca_window <- tile_context_chr_split(ca_split, "CA", s, depth, mincov, ca_cov_bases)
      metrics_list[[length(metrics_list)+1L]] <- compute_com_metrics(cg_window, ca_window, s, file_index$age[i], depth, mincov)
    }
  }
  metrics_dt <- rbindlist(metrics_list, use.names=TRUE, fill=TRUE)
  fwrite(metrics_dt, metric_cache, sep="\t")
}

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_end.txt"))
log_msg("02 complete")
log_msg(paste("Metric table:", metric_cache))
