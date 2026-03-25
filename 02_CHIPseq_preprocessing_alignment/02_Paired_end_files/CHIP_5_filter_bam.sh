#!/bin/bash

# 1. Environment Setup
conda activate methylation

# 2. Path Configuration
IN_DIR="./aln_dedup"
OUT_DIR="./processed_pe"
LOG_DIR="$OUT_DIR/logs"
THREADS=9

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

echo "--------------------------------------------------------"
echo "Starting Paired-End BAM Filtering Pipeline"
echo "--------------------------------------------------------"

# 3. Processing Loop
for bam_file in "$IN_DIR"/*_deduplicated.bam; do
    
    [ -e "$bam_file" ] || continue
    bname=$(basename "$bam_file" _deduplicated.bam)

    echo "========================================================"
    echo "Processing Sample: $bname"
    date
    echo "========================================================"

    # Step 1: Filter by MAPQ, proper pair, and alignment status
    # -q 10: MAPQ >= 10
    # -f 2: Keep only proper pairs
    # -F 268: Discard unmapped (4), mate unmapped (8), and not primary alignment (256)
    echo "Status: Filtering by MAPQ and Pair Status..."
    samtools view -@ "$THREADS" -q 10 -f 2 -F 268 -b -h "$bam_file" > "$OUT_DIR/${bname}_temp.bam"

    # Step 2: Index temporary BAM for chromosome filtering
    samtools index -@ "$THREADS" "$OUT_DIR/${bname}_temp.bam"

    # Step 3: Strict Chromosome Filtering
    # Regex targets standard mammalian chromosomes (chr1-22, X, Y, MT) to avoid scaffold artifacts
    echo "Status: Filtering out unassigned contigs/scaffolds..."
    samtools idxstats "$OUT_DIR/${bname}_temp.bam" | cut -f1 | grep -E '^chr([0-9]+|[XY]|MT|M)$' | \
        xargs samtools view -@ "$THREADS" -b -h -o "$OUT_DIR/${bname}_final.bam" "$OUT_DIR/${bname}_temp.bam"

    # Step 4: Index final BAM
    samtools index -@ "$THREADS" "$OUT_DIR/${bname}_final.bam"

    # Step 5: Clean up temporary files
    rm "$OUT_DIR/${bname}_temp.bam"
    rm "$OUT_DIR/${bname}_temp.bam.bai"

    # Step 6: Final Quality Control Stats
    echo "Status: Generating final alignment statistics..."
    samtools idxstats "$OUT_DIR/${bname}_final.bam" > "$LOG_DIR/${bname}_final_idxstats.txt"
    samtools flagstat -@ "$THREADS" "$OUT_DIR/${bname}_final.bam" > "$LOG_DIR/${bname}_final_flagstat.txt"
    samtools stats -@ "$THREADS" "$OUT_DIR/${bname}_final.bam" > "$LOG_DIR/${bname}_final_stats.txt"

    echo "Status: Finished $bname"
    echo "--------------------------------------------------------"

done

echo "Paired-End Filtering Completed. Clean BAM files are in $OUT_DIR"
