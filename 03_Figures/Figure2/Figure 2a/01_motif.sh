#!/bin/bash

# 1. Environment Setup
conda activate methylation

# 2. Output Directories
OUT_POS="homer_results_pos"
OUT_NEG="homer_results_neg"

mkdir -p "$OUT_POS"
mkdir -p "$OUT_NEG"

# 3. Motif Analysis for Positive Regions (e.g., Hypermethylated)
echo "Running motif analysis for COM_pos_motifs..."
findMotifsGenome.pl COM_pos_motifs hg38 "$OUT_POS" \
    -size given \
    -mask \
    -p 8 \
    -bits

# 4. Motif Analysis for Negative Regions (e.g., Hypomethylated)
echo "Running motif analysis for COM_neg_motifs..."
findMotifsGenome.pl COM_neg_motifs hg38 "$OUT_NEG" \
    -size given \
    -mask \
    -p 8 \
    -bits

echo "Motif enrichment analysis completed."
