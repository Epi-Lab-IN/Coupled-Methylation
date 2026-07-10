#!/bin/bash

# 1. Initialize Conda for shell script compatibility
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate methylation

# 2. Create output directories to keep the workspace clean
mkdir -p fastqc_pre_results

# 3. Perform Quality Control
echo "Running FastQC..."
fastqc *.fastq.gz -o fastqc_pre_results/ 

