#!/bin/bash

# 1. Environment Setup
conda activate methylation

# 2. Path Configuration
BAM_DIR=~/Escritorio/genome/aln
OUTDIR=~/Escritorio/results/aln_dedup
LOG_DIR="$OUTDIR/logs"  # Keeps the main folder clean from text files
THREADS=9

# Create output and log directories
mkdir -p "$OUTDIR"
mkdir -p "$LOG_DIR"

# 3. Processing Loop
# Iterate directly over all sorted BAM files in the input directory
for bam_file in "$BAM_DIR"/*_sorted.bam; do

    # Safety check: skip if no files match the pattern
    [ -e "$bam_file" ] || continue

    # Extract base name (removes path and _sorted.bam suffix)
    bname=$(basename "$bam_file" _sorted.bam)

    echo "========================================================"
    echo "STARTING DEDUPLICATION: $bname"
    date
    echo "========================================================"

    # Step 1: Mark and Remove Duplicates using Picard
    # REMOVE_DUPLICATES=true physically removes the duplicates instead of just flagging them.
    # This is standard for ChIP-seq peak calling to reduce false positive peaks.
    echo "Status: Running Picard MarkDuplicates..."
    picard MarkDuplicates \
        VERBOSITY=ERROR \
        REMOVE_DUPLICATES=true \
        ASSUME_SORT_ORDER=coordinate \
        MAX_RECORDS_IN_RAM=10000000 \
        INPUT="$bam_file" \
        METRICS_FILE="$LOG_DIR/${bname}_markdup_metrics.txt" \
        OUTPUT="$OUTDIR/${bname}_deduplicated.bam"

    # Step 2: Index the newly deduplicated BAM file
    echo "Status: Indexing deduplicated BAM..."
    samtools index -@ "$THREADS" "$OUTDIR/${bname}_deduplicated.bam"

    # Step 3: Generate Alignment Statistics
    # These metrics are crucial for the supplementary material (Mapping Rates, Duplicate Rates)
    echo "Status: Generating alignment statistics..."

    # idxstats: reports alignments per chromosome (useful to check mtDNA contamination)
    samtools idxstats "$OUTDIR/${bname}_deduplicated.bam" \
        > "$LOG_DIR/${bname}_idxstats.txt"

    # flagstat: reports basic mapping statistics
    samtools flagstat -@ "$THREADS" "$OUTDIR/${bname}_deduplicated.bam" \
        > "$LOG_DIR/${bname}_flagstat.txt"

    # stats: comprehensive statistics
    samtools stats -@ "$THREADS" "$OUTDIR/${bname}_deduplicated.bam" \
        > "$LOG_DIR/${bname}_stats.txt"

    echo "Status: Finished processing $bname"
    date
    echo "--------------------------------------------------------"

done

echo "Workflow completed. Deduplicated files are in $OUTDIR"
echo "QC metrics are safely stored in $LOG_DIR"
