# Extract average insert sizes from bam files
# using samtools stats, and save to a tsv file
# with accompaning sample names.

#!/bin/bash

# Directory containing BAM files
input_directory=$1

# Output TSV file
output_tsv=$2

# Write header to output TSV file
echo -e "Sample\tInsert_Size" > "$output_tsv"

# Iterate through BAM files in the directory
for bam_file in "$input_directory"/*.bam; do
    # Extract sample name from BAM file name
    sample_name=$(basename "$bam_file" .bam)

    # Extract insert sizes using samtools stats and parse average insert size
    avg_insert_line=$(samtools stats "$bam_file" | grep 'insert size average')
    insert_size=$(echo "$avg_insert_line" | awk '{print $5}')

    # Write sample name and insert size to output TSV file
    echo -e "$sample_name\t$insert_size" >> "$output_tsv"
done