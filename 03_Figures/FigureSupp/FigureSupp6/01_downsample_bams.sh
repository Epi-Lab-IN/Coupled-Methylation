#!/usr/bin/env bash
# 01_downsample_bams.sh
#
# PURPOSE
#   Sequentially downsample coordinate-sorted BAMs and run Bismark methylation extraction
#
#   This script writes a portable, consistent final file layout:
#     <EXTERNAL_ROOT>/{sample}/CX_processed/Data/
#       {sample}_{depth}x_unsorted.CX_report.txt_CG.gz
#       {sample}_{depth}x_unsorted.CX_report.txt_CA.gz
#       {sample}_{depth}x_unsorted.CX_report.txt_CC.gz
#
# SAFETY FEATURES FROM previous versions:
#   - target-specific work directories;
#   - free-space guard before starting/continuing each target;
#   - gzip validation of CX_report before filtering;
#   - gzip validation of final CG/CA/CC outputs;
#   - skip completed target only if final CG/CA/CC outputs are valid;
#   - restart incomplete target by removing target-specific work directory;
#   - do not mark completion when Bismark writes corrupt/truncated gzip files.
#
# IMPORTANT
#   samtools coverage is run without a thread option because some installed
#   samtools versions do not support multithreading for this subcommand.
# 	
# samtools and bismark must be installed within the environment

set -euo pipefail

############################
## User configuration
############################

# Edit these three paths for the current machine before running.
# BAM_ROOT contains <sample>_sorted.bam input files.
# EXTERNAL_ROOT is the portable 01-results directory written by this script.
# GENOME_FOLDER read by Bismark, contains the .fasta genome and a bisulfite-prepared genome subdirectory 
BAM_ROOT="/path/to/input_sorted_bam_directory"
EXTERNAL_ROOT="/path/to/01_results"
GENOME_FOLDER="/path/to/bismark_genome_folder/hg38_noXYM"

# Samples to process. Edit as needed.
SAMPLES=(M1 M2 M3 M4 M5 M6 M7)

# Target nominal mean sequencing depths.
TARGET_DEPTHS=(1 5 10 15 20 25)

# Mean depths from samtools-coverage-table.txt.
declare -A MEAN_DEPTH=(
  [M1]=26.1
  [M2]=31.5
  [M3]=26.5
  [M4]=26.4
  [M5]=21.2
  [M6]=25.2
  [M7]=21.9
)

SAMTOOLS_THREADS=8
BISMARK_MULTICORE=2
IGNORE_R2=4
SEED_BASE=202606

# External disk guard. Bismark --cytosine_report can require very large temporary
# files. Keep this high enough to prevent truncation/crashes.
MIN_EXTERNAL_FREE_GB=500

RUN_COVERAGE_QC=true
COVERAGE_AUTOSOMES_ONLY=true

# Restart behavior.
# If target work dir exists but final external CG/CA/CC files are missing/corrupt,
# remove target-specific work dir and retry.
FORCE_RESTART_INCOMPLETE_TARGET=true

# Skip only when final context outputs exist and pass gzip validation.
SKIP_COMPLETED=true

# Cleanup behavior on external drive.
# Coordinate BAM can be deleted after name-sort and coverage QC.
DELETE_COORD_BAM_AFTER_NAMESORT=true
# Name-sorted BAM can be deleted after methylation extraction succeeds.
DELETE_NAME_BAM_AFTER_EXTRACTION=true
# Full CX report is huge; delete after final CG/CA/CC have been validated.
DELETE_CX_REPORT_AFTER_FILTER=true
# Delete Bismark bedGraph/cov intermediates after final files are valid.
DELETE_BISMARK_LARGE_INTERMEDIATES=true
# Keep M-bias/splitting reports by copying them to qc_reports before cleaning.
COPY_BISMARK_QC_REPORTS=true

############################
## Optional conda activation
############################
# source "$HOME/miniconda3/etc/profile.d/conda.sh"
# conda activate bismark_env

############################
## Helpers
############################

log_msg() {
  local msg="$1"
  echo "$(date '+%Y-%m-%d %H:%M:%S') - ${msg}" | tee -a "$LOG_FILE"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: required command not found: $1" >&2; exit 1; }
}

free_gb_for_path() {
  local path="$1"
  mkdir -p "$path"
  df -BG "$path" | awk 'NR==2 {gsub("G","",$4); print $4}'
}

check_external_space() {
  local free_gb
  free_gb=$(free_gb_for_path "$EXTERNAL_ROOT")
  log_msg "External free space at ${EXTERNAL_ROOT}: ${free_gb} GB"
  if [[ "$free_gb" -lt "$MIN_EXTERNAL_FREE_GB" ]]; then
    log_msg "ERROR: external free space ${free_gb} GB is below MIN_EXTERNAL_FREE_GB=${MIN_EXTERNAL_FREE_GB}. Stop before target to avoid truncated gzip files."
    exit 1
  fi
}

fraction_for_target() {
  awk -v t="$1" -v m="$2" 'BEGIN { if (m <= 0) print "NA"; else printf "%.10f", t/m }'
}

leq_float() {
  awk -v a="$1" -v b="$2" 'BEGIN { exit !(a <= b) }'
}

samtools_subsample_arg() {
  local seed="$1"
  local frac="$2"
  echo "${seed}${frac#0}"
}

gzip_valid_nonempty() {
  local f="$1"
  [[ -s "$f" ]] && gzip -t "$f" >/dev/null 2>&1
}

final_contexts_valid() {
  local cg="$1" ca="$2" cc="$3"
  gzip_valid_nonempty "$cg" && gzip_valid_nonempty "$ca" && gzip_valid_nonempty "$cc"
}

write_global_coverage_header_if_needed() {
  local global_file="${EXTERNAL_ROOT}/all_samples_downsampled_bam_coverage_summary.tsv"
  if [[ ! -s "$global_file" ]]; then
    echo -e "sample\ttarget_depth_x\tconfigured_original_mean_depth_x\tsubsample_fraction\tsamtools_seed\tachieved_autosomal_mean_depth_x\tachieved_autosomal_coverage_pct\tautosomal_length_bp\tautosomal_covbases\tautosomal_numreads_sum\tcoverage_tsv" > "$global_file"
  fi
}

run_samtools_coverage_qc() {
  local bam="$1" sample="$2" target="$3" original_mean="$4" fraction="$5" seed="$6" cov_dir="$7"
  mkdir -p "$cov_dir"
  local depth_pad cov_tsv summary_tsv sample_all_tsv global_tsv
  depth_pad=$(printf "%02d" "$target")
  cov_tsv="${cov_dir}/${sample}_${depth_pad}x_samtools_coverage.tsv"
  summary_tsv="${cov_dir}/${sample}_${depth_pad}x_samtools_coverage.autosomes_summary.tsv"
  sample_all_tsv="${cov_dir}/${sample}_all_depths_samtools_coverage_summary.tsv"
  global_tsv="${EXTERNAL_ROOT}/all_samples_downsampled_bam_coverage_summary.tsv"

  if [[ -s "$summary_tsv" ]]; then
    log_msg "Coverage QC exists for ${sample} ${target}x: ${summary_tsv}"
    return 0
  fi

  log_msg "Indexing coordinate BAM for coverage QC: $bam"
  samtools index -@ "$SAMTOOLS_THREADS" "$bam"

  log_msg "Running samtools coverage for ${sample} ${target}x (single-threaded)"
  samtools coverage "$bam" -o "$cov_tsv"

  awk -v sample="$sample" -v target="$target" -v original_mean="$original_mean" -v fraction="$fraction" -v seed="$seed" -v cov_tsv="$cov_tsv" -v autosomes_only="$COVERAGE_AUTOSOMES_ONLY" '
    BEGIN { OFS="\t" }
    NR == 1 { next }
    {
      chr=$1; keep=1
      if (autosomes_only == "true") { c=chr; sub(/^chr/, "", c); keep=(c ~ /^[0-9]+$/ && c>=1 && c<=22) }
      if (keep) {
        len=$3-$2+1; numreads=$4; covbases=$5; meandepth=$7
        sum_len += len; sum_covbases += covbases; sum_depth_bases += meandepth*len; sum_numreads += numreads
      }
    }
    END {
      achieved_depth=(sum_len>0 ? sum_depth_bases/sum_len : "NA")
      achieved_cov_pct=(sum_len>0 ? 100*sum_covbases/sum_len : "NA")
      print "sample","target_depth_x","configured_original_mean_depth_x","subsample_fraction","samtools_seed","achieved_autosomal_mean_depth_x","achieved_autosomal_coverage_pct","autosomal_length_bp","autosomal_covbases","autosomal_numreads_sum","coverage_tsv"
      print sample,target,original_mean,fraction,seed,achieved_depth,achieved_cov_pct,sum_len,sum_covbases,sum_numreads,cov_tsv
    }' "$cov_tsv" > "$summary_tsv"

  if [[ ! -s "$sample_all_tsv" ]]; then head -n 1 "$summary_tsv" > "$sample_all_tsv"; fi
  tail -n 1 "$summary_tsv" >> "$sample_all_tsv"
  write_global_coverage_header_if_needed
  tail -n 1 "$summary_tsv" >> "$global_tsv"
  log_msg "Coverage QC summary: $(tail -n 1 "$summary_tsv")"
  rm -f "${bam}.bai" "${bam%.bam}.bai"
}

# splits cx_report into separate files containing only CG, CA, or CC dinucleotides
filter_cx_report() {
  local cx_report="$1"
  local out_prefix="$2"
  local final_dir="$3"

  mkdir -p "$final_dir"
  local cg_out="${final_dir}/${out_prefix}_CG.gz"
  local ca_out="${final_dir}/${out_prefix}_CA.gz"
  local cc_out="${final_dir}/${out_prefix}_CC.gz"

  log_msg "Filtering CX report into CG/CA/CC contexts: $(basename "$cx_report")"
  zcat "$cx_report" | awk 'BEGIN{OFS="\t"} ($6=="CG" || $7 ~ /^CG/) {print}' | gzip -c > "$cg_out"
  zcat "$cx_report" | awk 'BEGIN{OFS="\t"} ($7 ~ /^CA/) {print}' | gzip -c > "$ca_out"
  zcat "$cx_report" | awk 'BEGIN{OFS="\t"} ($7 ~ /^CC/) {print}' | gzip -c > "$cc_out"

  gzip -t "$cg_out"; gzip -t "$ca_out"; gzip -t "$cc_out"
  local n_cg n_ca n_cc
  n_cg=$(zcat "$cg_out" | wc -l); n_ca=$(zcat "$ca_out" | wc -l); n_cc=$(zcat "$cc_out" | wc -l)
  log_msg "Context rows: CG=${n_cg}; CA=${n_ca}; CC=${n_cc}"
  if [[ "$n_cg" -eq 0 || "$n_ca" -eq 0 || "$n_cc" -eq 0 ]]; then
    log_msg "ERROR: one or more context outputs has zero rows for ${out_prefix}"
    exit 1
  fi
}

find_cx_report() {
  find "$1" -maxdepth 1 -type f -name "*$2*CX_report.txt.gz" -size +0c | sort | tail -n 1
}

copy_qc_reports() {
  local extract_dir="$1" qc_dir="$2" bam_base="$3"
  [[ "$COPY_BISMARK_QC_REPORTS" == true ]] || return 0
  mkdir -p "$qc_dir"
  find "$extract_dir" -maxdepth 1 -type f \( \
    -name "*${bam_base}*M-bias.txt" -o \
    -name "*${bam_base}*splitting_report.txt" -o \
    -name "*${bam_base}*cytosine_context_summary.txt" \
  \) -print -exec cp -f {} "$qc_dir" \; | tee -a "$LOG_FILE" || true
}

cleanup_large_intermediates() {
  local extract_dir="$1" bam_base="$2"
  [[ "$DELETE_BISMARK_LARGE_INTERMEDIATES" == true ]] || return 0
  find "$extract_dir" -maxdepth 1 -type f \( \
    -name "*${bam_base}*.bedGraph.gz" -o \
    -name "*${bam_base}*.bismark.cov.gz" -o \
    -name "*${bam_base}*.coverage2cytosine.*" \
  \) -print -delete | tee -a "$LOG_FILE" || true
}

############################
## Check commands
############################
require_cmd samtools
require_cmd bismark_methylation_extractor
require_cmd awk
require_cmd gzip
require_cmd zcat

############################
## Main loop
############################
mkdir -p "$EXTERNAL_ROOT"

for SAMPLE in "${SAMPLES[@]}"; do
  INPUT_BAM="${BAM_ROOT}/${SAMPLE}_sorted.bam"
  SAMPLE_DIR="${EXTERNAL_ROOT}/${SAMPLE}"
  BAM_DOWN_DIR="${SAMPLE_DIR}/bam_downsampled"
  BAM_NAME_DIR="${SAMPLE_DIR}/bam_name_sorted"
  EXTRACT_ROOT="${SAMPLE_DIR}/meth_extract"
  FINAL_DIR="${SAMPLE_DIR}/CX_processed/Data"
  COVERAGE_DIR="${SAMPLE_DIR}/coverage_qc"
  QC_REPORT_DIR="${SAMPLE_DIR}/bismark_qc_reports"
  LOG_DIR="${SAMPLE_DIR}/logs"
  mkdir -p "$BAM_DOWN_DIR" "$BAM_NAME_DIR" "$EXTRACT_ROOT" "$FINAL_DIR" "$COVERAGE_DIR" "$QC_REPORT_DIR" "$LOG_DIR"
  LOG_FILE="${LOG_DIR}/${SAMPLE}_subsample_extract_filter_externalWork_v05.log"
  touch "$LOG_FILE"

  log_msg "===== Starting sample ${SAMPLE} ====="
  log_msg "Input BAM: ${INPUT_BAM}"
  check_external_space

  if [[ ! -s "$INPUT_BAM" ]]; then log_msg "ERROR: input BAM missing or empty: ${INPUT_BAM}"; continue; fi
  MEAN="${MEAN_DEPTH[$SAMPLE]:-}"
  if [[ -z "$MEAN" ]]; then log_msg "ERROR: no mean depth configured for ${SAMPLE}; skipping"; continue; fi

  if [[ ! -s "${INPUT_BAM}.bai" && ! -s "${INPUT_BAM%.bam}.bai" ]]; then
    log_msg "Indexing input BAM"
    samtools index -@ "$SAMTOOLS_THREADS" "$INPUT_BAM" || log_msg "WARNING: samtools index failed for input BAM"
  fi

  for TARGET in "${TARGET_DEPTHS[@]}"; do
    check_external_space
    if ! leq_float "$TARGET" "$MEAN"; then log_msg "SKIP ${SAMPLE} ${TARGET}x: target depth exceeds mean (${MEAN}x)"; continue; fi

    DEPTH_PAD=$(printf "%02d" "$TARGET")
    BASE="${SAMPLE}_${DEPTH_PAD}x_unsorted"
    FINAL_CG="${FINAL_DIR}/${BASE}.CX_report.txt_CG.gz"
    FINAL_CA="${FINAL_DIR}/${BASE}.CX_report.txt_CA.gz"
    FINAL_CC="${FINAL_DIR}/${BASE}.CX_report.txt_CC.gz"

    if [[ "$SKIP_COMPLETED" == true ]] && final_contexts_valid "$FINAL_CG" "$FINAL_CA" "$FINAL_CC"; then
      log_msg "SKIP completed ${SAMPLE} ${TARGET}x: valid final CG/CA/CC files exist"
      continue
    fi

    TARGET_WORK_DIR="${SAMPLE_DIR}/work_${BASE}"
    EXTRACT_DIR="${TARGET_WORK_DIR}/meth_extract"
    mkdir -p "$TARGET_WORK_DIR" "$EXTRACT_DIR"

    if [[ -d "$TARGET_WORK_DIR" && "$FORCE_RESTART_INCOMPLETE_TARGET" == true && ! -s "$FINAL_CG" ]]; then
      log_msg "Removing previous incomplete target work directory: $TARGET_WORK_DIR"
      rm -rf "$TARGET_WORK_DIR"
      mkdir -p "$TARGET_WORK_DIR" "$EXTRACT_DIR"
    fi

    FRAC=$(fraction_for_target "$TARGET" "$MEAN")
    SEED=$((SEED_BASE + TARGET))
    SUBARG=$(samtools_subsample_arg "$SEED" "$FRAC")
    COORD_BAM="${TARGET_WORK_DIR}/${SAMPLE}_${DEPTH_PAD}x_coord.bam"
    NAME_BAM="${TARGET_WORK_DIR}/${BASE}.bam"

    log_msg "----- ${SAMPLE} target ${TARGET}x: fraction=${FRAC}, samtools -s ${SUBARG} -----"
    log_msg "External target work directory: $TARGET_WORK_DIR"

    log_msg "Subsampling coordinate-sorted BAM on external drive"
    samtools view -@ "$SAMTOOLS_THREADS" -s "$SUBARG" -b -o "$COORD_BAM" "$INPUT_BAM"

    if [[ "$RUN_COVERAGE_QC" == true ]]; then
      run_samtools_coverage_qc "$COORD_BAM" "$SAMPLE" "$TARGET" "$MEAN" "$FRAC" "$SEED" "$COVERAGE_DIR"
    fi

    log_msg "Name-sorting BAM for paired-end methylation extraction"
    samtools sort -n -@ "$SAMTOOLS_THREADS" -o "$NAME_BAM" "$COORD_BAM"
    if [[ "$DELETE_COORD_BAM_AFTER_NAMESORT" == true ]]; then
      rm -f "$COORD_BAM" "${COORD_BAM}.bai" "${COORD_BAM%.bam}.bai"
    fi

    check_external_space
    log_msg "Running bismark_methylation_extractor with external target output_dir"
    bismark_methylation_extractor \
      --paired-end --no_overlap --gzip --bedGraph --cytosine_report -CX \
      --genome_folder "$GENOME_FOLDER" \
      --multicore "$BISMARK_MULTICORE" \
      --output_dir "$EXTRACT_DIR" \
      --ignore_r2 "$IGNORE_R2" \
      "$NAME_BAM" 2>&1 | tee -a "$LOG_FILE"

    BAM_BASE=$(basename "$NAME_BAM" .bam)
    CX_REPORT=$(find_cx_report "$EXTRACT_DIR" "$BAM_BASE")
    if [[ -z "$CX_REPORT" || ! -s "$CX_REPORT" ]]; then
      log_msg "ERROR: no non-empty CX_report.txt.gz found for ${BAM_BASE} in ${EXTRACT_DIR}"
      exit 1
    fi
    if ! gzip_valid_nonempty "$CX_REPORT"; then
      log_msg "ERROR: CX_report is corrupt or incomplete: $CX_REPORT"
      exit 1
    fi

    filter_cx_report "$CX_REPORT" "${BASE}.CX_report.txt" "$FINAL_DIR"
    if ! final_contexts_valid "$FINAL_CG" "$FINAL_CA" "$FINAL_CC"; then
      log_msg "ERROR: final context files failed gzip validation"
      exit 1
    fi

    copy_qc_reports "$EXTRACT_DIR" "$QC_REPORT_DIR" "$BAM_BASE"
    cleanup_large_intermediates "$EXTRACT_DIR" "$BAM_BASE"

    if [[ "$DELETE_NAME_BAM_AFTER_EXTRACTION" == true ]]; then
      rm -f "$NAME_BAM"
    fi
    if [[ "$DELETE_CX_REPORT_AFTER_FILTER" == true ]]; then
      rm -f "$CX_REPORT"
    fi

    # Remove target work directory if empty after cleanup.
    rmdir "$EXTRACT_DIR" 2>/dev/null || true
    rmdir "$TARGET_WORK_DIR" 2>/dev/null || true

    check_external_space
    log_msg "DONE ${SAMPLE} ${TARGET}x"
  done

  log_msg "===== Finished sample ${SAMPLE} ====="
done

echo "All requested samples/depths processed."
