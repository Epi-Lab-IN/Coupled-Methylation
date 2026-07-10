#!/bin/bash/

# 1. Environment Setup
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate methylation

# 2. Configuration
# Path to the reference FASTA file
GENOME_FASTA="/path/to/genome.fa"

# Path and Base Name for the output index
# Example: ./index/hg38
OUTPUT_INDEX_BASE="/path/to/folder/proyect/hg38_index"


# 3. Run bowtie2-build
bowtie2-build --threads 10 "$GENOME_FASTA" "$OUTPUT_INDEX_BASE"
