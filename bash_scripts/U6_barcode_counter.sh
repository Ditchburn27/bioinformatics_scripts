#! /bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <barcodes.csv> <reads.fastq.gz>"
    exit 1
fi

# Input files
csv_file="$1"
fastq_file="$2"

# Extract library_id from the fastq file name (assumes it's the first part before the first underscore)
library_id=$(basename "$fastq_file" | cut -d_ -f1)

# Output file name using the library_id
output_file="${library_id}_barcode_counts.csv"

# Calculate the total number of reads
total_reads=$(zcat "$fastq_file" | awk 'NR%4==2' | wc -l)

# Create a temporary file to store intermediate results
tmp_file=$(mktemp)

# Process each barcode and write results to temporary file
awk -F, -v total_reads="$total_reads" '
NR>1 {
    cmd="zcat \"'$fastq_file'\" | egrep \""$1"\" | awk \"NR%4==2\" | wc -l";
    cmd | getline count; close(cmd);
    sum+=count;
    printf("%s,%s,%s,%d,%.2f\n", $1, $2, $3, count, (count/total_reads)*100) >> "tmp_file"
}
END {total_barcode_counts=sum}
' "$csv_file"

# Calculate total barcode reads separately
total_barcode_counts=$(awk -F, '{sum+=$4} END {print sum}' tmp_file)

# Add header and initial rows to the output file
{
    echo "Total Reads,Barcode Reads"
    echo "$total_reads,$total_barcode_counts"
    echo ""
    echo "Barcode,Name,TF,Counts,Barcode proportion of total reads (%),Barcode proportion of barcode reads (%)"
} > "$output_file"

# Append data with proportions to the output file
awk -F, -v total_barcode_counts="$total_barcode_counts" '
{
    
    print $0 "," sprintf("%.2f", ($4/total_barcode_counts)*100)
    
}
' tmp_file >> "$output_file"

# Cleanup temporary file
rm tmp_file