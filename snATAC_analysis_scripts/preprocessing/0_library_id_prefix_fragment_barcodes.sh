#!/bin/bash

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 input_file.tsv.gz"
    exit 1
fi

# Process input file
gzip -dc "$1" | awk -F"\t" -v OFS="\t" '
    BEGIN {flag=0} 
    { 
        if ($1 ~ /^# id=/) {
            parent=substr($1, 7)
            print $0  # Preserve the library ID header line
            next
        }
        if (flag==0) {
            if ($1 ~ /^chr/) {
                flag=1
            }
        } 
        if (flag==1 && $1 ~ /^chr/) {
            $4 = "R" parent "+" $4  # Prepend library ID to the barcode
        }
        print
    }' > output_file.tsv

echo "Processing complete. Output saved to output_file.tsv"