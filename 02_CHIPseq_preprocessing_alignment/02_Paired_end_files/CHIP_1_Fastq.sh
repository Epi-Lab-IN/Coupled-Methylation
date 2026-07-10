#!/bin/bash/
# 1. Environment Setup
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate trim

# 2. Directory Configuration
INPUT_DIR="./raw_data_geo"
OUTPUT_DIR="./qc_results"

mkdir -p "$OUTPUT_DIR"

# 3. Processing Loop
# This handles paired-end files 
for R1 in "$INPUT_DIR"/*_1.fastq.gz; do
    
    [ -e "$R1" ] || continue
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    SAMPLE_NAME=$(basename "$R1" _1.fastq.gz)

    echo "Processing Sample: $SAMPLE_NAME"

    # Running Trim Galore with the specific parameters mentioned:
    # -l 20: Length filtering
    # -p / --paired: Paired-end mode
    # --fastqc: Run FastQC automatically after trimming
    
    trim_galore --paired \
        --fastqc \
        --gzip \
        --length 20 \
        --output_dir "$OUTPUT_DIR" \
        "$R1" "$R2"

done

echo "--------------------------------------------------------"
echo "QC and Preprocessing completed. Results in $OUTPUT_DIR"
