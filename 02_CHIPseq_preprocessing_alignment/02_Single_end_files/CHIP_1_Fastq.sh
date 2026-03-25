#!/bin/bash/
# 1. Environment Setup
conda activate trim

# 2. Directory Configuration
INPUT_DIR="./raw_data_geo"
OUTPUT_DIR="./qc_results"

mkdir -p "$OUTPUT_DIR"

# 3. Processing Loop
# This handles single-end files 
for FASTQ in "$INPUT_DIR"/*.fastq.gz; do
    
    [ -e "$FASTQ" ] || continue
    SAMPLE_NAME=$(basename "$FASTQ" .fastq.gz)

    echo "Processing Sample: $SAMPLE_NAME"

    # Trim Galore parameters for Single-End:
    # --fastqc: Runs FastQC on the fly
    # --length 20: Filters out reads shorter than 20bp (Ref: -l 20)
    
    trim_galore \
        --fastqc \
        --gzip \
        --length 20 \
        --output_dir "$OUTPUT_DIR" \
        "$FASTQ"

done

echo "--------------------------------------------------------"
echo "Single-End QC completed. Results saved in $OUTPUT_DIR"
