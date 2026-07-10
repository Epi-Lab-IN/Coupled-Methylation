#!/bin/bash

# 1. Initialize Conda for shell script compatibility
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate trim

# 2. Create output directories to keep the workspace clean
mkdir -p trimmed_reads

# 3. Trim adapters and filter by quality (Paired-End)
echo "Starting Trim Galore..."
trim_galore --paired -o trimmed_reads *_R1*.fastq.gz *_R2*.fastq.gz
