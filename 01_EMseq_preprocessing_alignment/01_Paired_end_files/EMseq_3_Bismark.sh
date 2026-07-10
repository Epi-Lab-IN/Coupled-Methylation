#!/bin/bash

# 1. Setup
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate methylation

GENOME_FOLDER="./genome"
THREADS=4 # CPUs per sample
mkdir -p alignment_results

# Prepare genome
echo "------------------------------------------"
echo "STEP 1: Preparing Genome"
echo "------------------------------------------"
bismark_genome_preparation --bowtie2 $GENOME_FOLDER

# 2. Loop through all R1 files in the trimmed_reads folder
for R1 in trimmed_reads/*_R1_val_1.fq.gz
do
    # Detect R2 by replacing R1 with R2 in the filename
    R2=${R1/_R1_val_1.fq.gz/_R2_val_2.fq.gz}
    
    # Extract sample name for logs (e.g., "Sample_A")
    SAMPLE_NAME=$(basename $R1 _R1_val_1.fq.gz)

    echo "=========================================="
    echo "PROCESSING SAMPLE: $SAMPLE_NAME"
    echo "=========================================="

    # 3. Alignment
    bismark --bowtie2 \
        --multicore $THREADS \
        --genome $GENOME_FOLDER \
        -1 $R1 -2 $R2 \
        --output_dir ./alignment_results/

    # 4. Deduplication
    deduplicate_bismark \
    	--bam \
        --paired \
        ./alignment_results/${SAMPLE_NAME}_R1_val_1_bismark_bt2_pe.bam

    # 5. Methylation Extraction (For paired-end)
    bismark_methylation_extractor \
        --paired \
        --no_overlap \
        --bedGraph \
        --cytosine_report\
        --CX\
        --gzip\
        ./alignment_results/${SAMPLE_NAME}_R1_val_1_bismark_bt2_pe.deduplicated.bam \
        --genome_folder $GENOME_FOLDER\
        -o ./alignment_results/

   
    echo "DONE: $SAMPLE_NAME"
done

echo "All samples processed."
