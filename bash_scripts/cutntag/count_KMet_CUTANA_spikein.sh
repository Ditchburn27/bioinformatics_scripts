#!/bin/bash
##########################################
# Written by Dr. Bryan Venters, EpiCypher Inc.
# Updated: Modified by [Leighton Ditchburn], [28.05.2025]
#
# Purpose: Use "grep -c" to count exact match to CUTANA spike-in nucleosome barcodes from unzipped paired-end (R1 & R2) fastqs. 
# This script has been modified to process all fastq files in a specified directory and output results to a CSV file.
##########################################

# Check if directory argument is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <directory_with_fastq_files> <output_csv_file>"
  exit 1
fi

# Input directory with fastq files and output CSV file
input_dir="$1"
output_csv="$2"

# List of barcodes
barcodes=("TTCGCGCGTAACGACGTACCGT" "CGCGATACGACCGCGTTACGCG" "CGACGTTAACGCGTTTCGTACG" "CGCGACTATCGCGCGTAACGCG" \
"CCGTACGTCGTGTCGAACGACG" "CGATACGCGTTGGTACGCGTAA" "TAGTTCGCGACACCGTTCGTCG" "TCGACGCGTAAACGGTACGTCG" \
"TTATCGCGTCGCGACGGACGTA" "CGATCGTACGATAGCGTACCGA" "CGCATATCGCGTCGTACGACCG" "ACGTTCGACCGCGGTCGTACGA" \
"ACGATTCGACGATCGTCGACGA" "CGATAGTCGCGTCGCACGATCG" "CGCCGATTACGTGTCGCGCGTA" "ATCGTACCGCGCGTATCGGTCG" \
"CGTTCGAACGTTCGTCGACGAT" "TCGCGATTACGATGTCGCGCGA" "ACGCGAATCGTCGACGCGTATA" "CGCGATATCACTCGACGCGATA" \
"CGCGAAATTCGTATACGCGTCG" "CGCGATCGGTATCGGTACGCGC" "GTGATATCGCGTTAACGTCGCG" "TATCGCGCGAAACGACCGTTCG" \
"CCGCGCGTAATGCGCGACGTTA" "CCGCGATACGACTCGTTCGTCG" "GTCGCGAACTATCGTCGATTCG" "CCGCGCGTATAGTCCGAGCGTA" \
"CGATACGCCGATCGATCGTCGG" "CCGCGCGATAAGACGCGTAACG" "CGATTCGACGGTCGCGACCGTA" "TTTCGACGCGTCGATTCGGCGA")

# Prepare the output CSV file with headers
echo -n "Barcode" > "$output_csv"
for file in "$input_dir"/*.fastq; do
  echo -n ",$(basename "$file")" >> "$output_csv"
done
echo "" >> "$output_csv"

# Iterate over each barcode and count occurrences in each file
for barcode in "${barcodes[@]}"; do
  echo -n "$barcode" >> "$output_csv"
  for file in "$input_dir"/*.fastq; do
    count=$(grep -c "$barcode" "$file")
    echo -n ",$count" >> "$output_csv"
  done
  echo "" >> "$output_csv"
done

echo "Barcode counts have been written to $output_csv."