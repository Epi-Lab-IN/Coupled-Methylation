#!/bin/bash

# 1. Environment Setup
conda activate epic2_env

# 2. Configuration & Paths
# Define standard parameters
BIN_SIZE=200
GAPS_ALLOWED=3
FDR_CUTOFF=0.05
MAPQ=10
GENOME=mm10

# Directories (Linked to your Paired-End output)
IN_DIR="./processed_pe"
OUT_DIR="./peaks_epic2_pe"
LOG_DIR="$OUT_DIR/logs"

# Ensure you specify exactly which BAM is your Input/Control
CONTROL_BAM="$IN_DIR/Input_final.bam" 
CHROM_SIZES="./reference_data/custom_genome.chromsizes"

# Create output directories
mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p $(dirname "$CHROM_SIZES")

echo "--------------------------------------------------------"
echo "STEP 1: Generating standard chrom.sizes from Control BAM"
echo "--------------------------------------------------------"

# Extract chromosome sizes directly from the BAM header using Perl
samtools view -H "$CONTROL_BAM" | perl -ne 'if(/^@SQ.*?SN:(\S+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > "$CHROM_SIZES"

echo "chrom.sizes successfully saved to $CHROM_SIZES"
echo "--------------------------------------------------------"
echo "STEP 2: Running epic2 Peak Calling (Paired-End mode)"
echo "--------------------------------------------------------"

# 3. Processing Loop for all Treatment BAMs
for TREATMENT_BAM in "$IN_DIR"/*_final.bam; do
    
    # Safety check: skip if file doesn't exist
    [ -e "$TREATMENT_BAM" ] || continue
    
    # Skip the control file so we don't run the Input against itself
    if [[ "$TREATMENT_BAM" == "$CONTROL_BAM" ]]; then
        continue
    fi

    # Extract base name (e.g., SampleA_final -> SampleA)
    bname=$(basename "$TREATMENT_BAM" _final.bam)

    echo "Processing Treatment: $bname vs Control"
    
    # Run epic2
    # --guess-bampe: CRUCIAL for Paired-End. Uses actual fragment lengths.
    epic2 \
        --treatment "$TREATMENT_BAM" \
        --control "$CONTROL_BAM" \
        --genome "$GENOME" \
        --keep-duplicates \
        --false-discovery-rate-cutoff "$FDR_CUTOFF" \
        --mapq "$MAPQ" \
        --chromsizes "$CHROM_SIZES" \
        --guess-bampe \
        --bin-size "$BIN_SIZE" \
        --gaps-allowed "$GAPS_ALLOWED" \
        --output "$OUT_DIR/${bname}_epic2_peaks.bed" \
        2> "$LOG_DIR/log_${bname}_epic2.txt"

    echo "Finished: $bname"

done

echo "--------------------------------------------------------"
echo "Peak calling completed. Broad BED files are available in $OUT_DIR"
