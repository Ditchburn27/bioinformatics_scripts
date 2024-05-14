#!/bin/bash

# Define the CSV file
csv_file="barcode_consensus_sequences.csv"

# Add CSV header if the file doesn't exist
if [ ! -f "$csv_file" ]; then
    echo "Library,Sequence,GC_Content" > "$csv_file"
fi

# Recursively search for .fa files in the results directory
find results -type f -name "*.fa" | while read -r fa_file; do
    # Extract library identifier from directory name
    library=$(basename "$(dirname "$fa_file")" | cut -d '_' -f 1)

    # Extract sequences and calculate GC content
    awk 'NR%4==2' "$fa_file" | cut -c 29-48 | while read -r sequence; do
        gc_content=$(printf "%.2f" "$(echo "scale=2; ($(grep -o '[GC]' <<< "$sequence" | wc -l) / ${#sequence}) * 100" | bc)")
        # Append results to CSV file
        echo "$library,$sequence,$gc_content" >> "$csv_file"
    done
done
