#!/bin/bash

# Check if dir_list.txt exists
if [ ! -f "dir_list.txt" ]; then
    echo "Error: dir_list.txt not found."
    exit 1
fi

# Concatenate all processed TSV files into one, avoiding duplicate headers
first_file=true
for directory in $(cat dir_list.txt); do
    file="$directory/outs/output_file.tsv"
    if [ -f "$file" ]; then
        if $first_file; then
            # Include header for the first file
            cat "$file" > temp_combined.tsv
            first_file=false
        else
            # Skip header for subsequent files
            grep -v '^#' "$file" >> temp_combined.tsv
        fi
    fi
done


mv temp_combined.tsv braindev_combined_fragments.tsv

echo "All fragments have been combined into braindev_combined_fragments.tsv"

## After running script remove extra contigs from fragments file:
#grep -v -f <(sed 's/^/^/' contig_list.txt) braindev_combined_fragments.tsv > filtered_combined_fragments.tsv

