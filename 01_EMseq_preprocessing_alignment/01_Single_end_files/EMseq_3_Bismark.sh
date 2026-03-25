#!/bin/bash

# 1. Setup
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
for R1 in trimmed_reads/*_trimmed.fq.gz
do
    
    # Extract sample name 
    SAMPLE_NAME=$(basename "$R1" _trimmed.fq.gz)

    echo "=========================================="
    echo "PROCESSING SAMPLE: $SAMPLE_NAME"
    echo "=========================================="

    # 3. Alignment
    bismark --bowtie2 \
        --multicore $THREADS \
        --genome $GENOME_FOLDER \
         "$R1" \
        --output_dir ./alignment_results/

    # 4. Deduplication
    deduplicate_bismark \
    	--bam \
        ./alignment_results/${SAMPLE_NAME}_trimmed_bismark_bt2.bam

    # 5. Methylation Extraction (For single-end)
    bismark_methylation_extractor \
        --bedGraph \
        --cytosine_report\
        --CX\
        --gzip\
        ./alignment_results/${SAMPLE_NAME}_trimmed_bismark_bt2.deduplicated.bam \
        --genome_folder $GENOME_FOLDER\
        -o ./alignment_results/

   
    echo "DONE: $SAMPLE_NAME"
done

echo "All samples processed."
