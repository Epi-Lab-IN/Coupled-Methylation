#!/bin/bash

# 1. Environment Setup
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate methylation

# 2. Path Configuration
# Specify the folder where your FASTQ files are
FASTQ_DIR="$HOME/Escritorio/raw_data"
GENOME_INDEX="$HOME/Escritorio/hg38_index/hg38" # Make sure this is the prefix!
OUTDIR="$HOME/Escritorio/results/aln"
THREADS=10

# Create output directory if it does not exist
mkdir -p "$OUTDIR"

echo "--------------------------------------------------------"
echo "Starting Paired-End Alignment Pipeline"
echo "--------------------------------------------------------"

# 3. Processing Loop
# We iterate ONLY over the R1 files to avoid processing pairs twice.
# Adjust the '*_R1.fastq.gz' pattern if your files end in '_1.fastq.gz'
for R1 in "$FASTQ_DIR"/*_R1.fastq.gz; do
    
    # Safety check: skip if no files match the pattern
    [ -e "$R1" ] || continue

    # Dynamically define the R2 file by replacing _R1 with _R2 in the string
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

    # Check if the R2 mate actually exists before aligning
    if [ ! -f "$R2" ]; then
        echo "ERROR: Mate file missing for $R1. Skipping..."
        continue
    fi

    # Define the base name for the sample dynamically (no hardcoded string cutting)
    finalname=$(basename "$R1" _R1.fastq.gz)

    echo "========================================================"
    echo "STARTING ALIGNMENT: $finalname"
    date
    echo "========================================================"

    # Step 1: Alignment with Bowtie2 (Paired-End Mode)
    # -1: Forward reads, -2: Reverse reads
    # -X 2000: Crucial for PE chromatin data to allow fragments up to 2kb
    # Results are piped directly to samtools to avoid large intermediate SAM files
    echo "Status: Running Bowtie2 Paired-End alignment..."
    bowtie2 -q --threads "$THREADS" -X 2000 --phred33 --very-sensitive-local \
        -x "$GENOME_INDEX" \
        -1 "$R1" -2 "$R2" \
        2> "$OUTDIR/log_${finalname}_aln.log" | \
        samtools view -@ "$THREADS" -bS - > "$OUTDIR/${finalname}_temp.bam"

    # Step 2: Sorting the BAM file
    echo "Status: Sorting BAM file..."
    samtools sort -@ "$THREADS" -o "$OUTDIR/${finalname}_aln_sorted.bam" "$OUTDIR/${finalname}_temp.bam"

    # Step 3: Indexing the BAM file
    echo "Status: Indexing BAM file..."
    samtools index -@ "$THREADS" "$OUTDIR/${finalname}_aln_sorted.bam"

    # Step 4: Cleanup
    # Remove the unsorted temporary BAM file to save space
    rm "$OUTDIR/${finalname}_temp.bam"

    echo "Status: Finished processing $finalname"
    date
    echo "--------------------------------------------------------"

done

echo "Workflow completed. All Paired-End samples aligned and indexed in $OUTDIR"
