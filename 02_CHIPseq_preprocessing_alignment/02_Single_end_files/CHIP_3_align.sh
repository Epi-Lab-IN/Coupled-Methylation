#!/bin/bash

# 1. Environment Setup
conda activate methylation

# 2. Path Configuration
# Define input FASTQ files, genome index prefix, and output directory
FASTQ_FILES=(~/Escritorio/*.fastq.gz)
GENOME_INDEX=~/Escritorio/hg38_index
OUTDIR=~/Escritorio/results/aln
THREADS=10

# Create output directory if it does not exist
mkdir -p "$OUTDIR"

# 3. Processing Loop
for fq in "${FASTQ_FILES[@]}"; do
    

    # Define the base name for the sample
    finalname=$(basename "$fq" .fastq.gz)

    echo "========================================================"
    echo "STARTING ALIGNMENT: $finalname"
    date
    echo "========================================================"

    # Step 1: Alignment with Bowtie2
    # --very-sensitive-local: high sensitivity while allowing soft-clipping
    # Results are piped directly to samtools to avoid large intermediate SAM files
    echo "Status: Running Bowtie2..."
    bowtie2 -q --threads "$THREADS" --phred33 --very-sensitive-local \
        -x "$GENOME_INDEX" \
        -U "$fq" \
        2> "$OUTDIR/log_${finalname}_aln.log" | \
        samtools view -@ "$THREADS" -bS - > "$OUTDIR/${finalname}_temp.bam"

    # Step 2: Sorting the BAM file
    # Required for downstream analysis and indexing
    echo "Status: Sorting BAM file..."
    samtools sort -@ "$THREADS" -o "$OUTDIR/${finalname}_aln_sorted.bam" "$OUTDIR/${finalname}_temp.bam"

    # Step 3: Indexing the BAM file
    # Generates .bai file for rapid coordinate-based access (IGV, deepTools)
    echo "Status: Indexing BAM file..."
    samtools index -@ "$THREADS" "$OUTDIR/${finalname}_aln_sorted.bam"

    # Step 4: Cleanup
    # Remove the unsorted temporary BAM file to save space
    rm "$OUTDIR/${finalname}_temp.bam"

    echo "Status: Finished processing $finalname"
    date
    echo "--------------------------------------------------------"

done

echo "Workflow completed. All samples aligned and indexed in $OUTDIR"
