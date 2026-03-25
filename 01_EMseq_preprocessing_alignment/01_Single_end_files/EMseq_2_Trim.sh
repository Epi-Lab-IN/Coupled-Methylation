#!/bin/bash

# 1. Initialize Conda for shell script compatibility
conda activate trim

# 2. Create output directories to keep the workspace clean
mkdir -p trimmed_reads


# 3 Trim adapters and filter by quality (Single-End)
echo "Starting Trim Galore..."
trim_galore -o trimmed_reads *.fastq.gz
