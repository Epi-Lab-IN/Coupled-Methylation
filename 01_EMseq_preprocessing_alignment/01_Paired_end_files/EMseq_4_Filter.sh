#!/bin/bash
# 1. Directory Configuration
INPUT_DIR="."
OUTPUT_DIR="./methylation_contexts"
# 2. Create output directory
mkdir -p "$OUTPUT_DIR"
echo "Starting context-specific extraction from CX reports..."
echo "Input Directory: $INPUT_DIR"
echo "Output Directory: $OUTPUT_DIR"
# 3. Processing Loop
for file in "$INPUT_DIR"/*.CX_report.txt.gz; do
    
    # Check if files exist
    [ -e "$file" ] || { echo "No .CX_report.txt.gz files found."; exit 1; }
    # Extract sample name (removes the long Bismark suffix)
    filename=$(basename "$file" .CX_report.txt.gz)
    echo "--------------------------------------------------------"
    echo "Processing Sample: $filename"
    # Core logic using pigz and awk
    # The CX report format:
    # col 1: chromosome | col 2: position | col 3: strand
    # col 4: count methylated | col 5: count unmethylated
    # col 6: context (CpG, etc) | col 7: trinucleotide context
    
    pigz -dc -p $(nproc) "$file" | awk -F'\t' -v out_dir="$OUTPUT_DIR" -v base="$filename" '
    BEGIN {
        OFS = "\t"
        cg_out  = "pigz -p 2 > " out_dir "/" base "_CG.gz"
        cc_out  = "pigz -p 2 > " out_dir "/" base "_CC.gz"
        ca_out  = "pigz -p 2 > " out_dir "/" base "_CA.gz"
    }
    {
        # 1. Filter: Skip non-informative positions (zero coverage)
        if ($4 == 0 && $5 == 0) next

        # 2. Filter: Skip CT context
        if ($7 ~ /^CT[ACGT]$/) next
        
        # 3. Redirect based on trinucleotide context (col 7)
        if ($7 ~ /^CG[ACGT]$/) {
            print | cg_out
        } else if ($7 ~ /^CC[ACGT]$/) {
            print | cc_out
        } else if ($7 ~ /^CA[ACGT]$/) {
            $6 = "CHH"
            print | ca_out
        }
    }
    END {
        close(cg_out)
        close(cc_out)
        close(ca_out)
    }
    '
    
    if [ $? -eq 0 ]; then
        echo "Done: $filename (CG, CC, and CA contexts exported)"
    else
        echo "Error processing: $filename"
    fi
done
echo "--------------------------------------------------------"
echo "Workflow completed. Output files stored in $OUTPUT_DIR"
