#!/bin/bash

sample=$1
barcode_type=$2
python_script=$3

UPS2_constant_region='GTCACTGTAAAACGACGGCCAGACTAGT'
RANDOM_BC_constant_region='CGTACGGCTTTAAGGCCGGTCCTAGCAA'
SERLOIN_constant_region='CGCCTCCCGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCGA'
BORG_constant_region='TGCAGCAGAGCCCGTAAGCCCGATAGCCC'

# Extract library ID from sample name
library_id=$(echo "$sample" | sed -E 's/.*(UBL[0-9]+)_.*/\1/')

# Create directory for the library
mkdir -p "${library_id}_bowtie_indexes"

if [ $barcode_type == "borg" ]; then
    constant_region=$BORG_constant_region
    reference_filename="/scratch/u6_barcode_fastqs/BORG_RNLS.fa"
    barcode_pos='29-48'
elif [ $barcode_type == "serloin" ]; then
    constant_region=$SERLOIN_constant_region
    reference_filename="/scratch/u6_barcode_fastqs/SERLOIN_RNLS.fa"
    barcode_pos='29-48'
elif [ $barcode_type == "random" ]; then
    constant_region=$RANDOM_BC_constant_region
    reference_filename="/scratch/u6_barcode_fastqs/U6_random_barcode_50_NT.fa"
    barcode_pos='29-78'
else
    echo "Invalid barcode type. Options are 'borg', 'serloin', or 'random'."
    exit 1
fi

# Extract barcodes
echo "Extracting barcodes with $barcode_type constant region..."
zcat $sample | egrep $UPS2_constant_region | egrep $constant_region | awk 'NR%4==2' | cut -c $barcode_pos > "${library_id}_output.txt"

# Generate reference fasta files with consensus barcodes
updated_reference_filename="${library_id}_${barcode_type}_reference.fa"
python "$python_script" "${library_id}_output.txt" "$reference_filename" "$barcode_type" "${library_id}"

# Index reference fasta file
bowtie2-build $updated_reference_filename "${library_id}_bowtie_indexes/${updated_reference_filename%.fa}" &> "${library_id}_bt2_index.log"

# Align paired-end reads to the indexed reference
bowtie2 -x "${library_id}_bowtie_indexes/${updated_reference_filename%.fa}" -1 "trimmed_fastqs/${sample/R2/R1}" -2 "trimmed_fastqs/${sample}" -S "${library_id}_alignment.sam" &> "${library_id}_bowtie.log"

# Convert SAM to BAM and sort
samtools view -b -o "${library_id}_alignment.bam" "${library_id}_alignment.sam"
samtools sort -o "${library_id}_alignment.sorted.bam" "${library_id}_alignment.bam"

# Index sorted BAM file
samtools index "${library_id}_alignment.sorted.bam"

# Generate alignment statistics summary
samtools stats "${library_id}_alignment.sorted.bam" > "${library_id}_alignment_stats.txt"

# Clean up intermediate files
rm -r "${library_id}_bowtie_indexes"
rm "${library_id}_alignment.sam" "${library_id}_alignment.bam" "${library_id}_output.txt"

# Make directory to store outputs
mkdir -p "${library_id}_results"
mv "${library_id}_alignment_stats.txt" "${library_id}_alignment.sorted.bam.bai" "${library_id}_alignment.sorted.bam" "${updated_reference_filename}" "${library_id}_bt2_index.log" "${library_id}_bowtie.log" "${library_id}_results/"

