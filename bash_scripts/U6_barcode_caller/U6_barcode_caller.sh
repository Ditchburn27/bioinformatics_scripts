#!/bin/bash

input_dir=$1
script_dir=$2
sampleID=$3

# change directory
cd $input_dir

# Trim fastqs in input_dir
bash ${script_dir}/fastp_trim_fastqs.sh $input_dir $sampleID

# BORG barcodes
parallel bash "${script_dir}/get_U6_barcodes.sh" {} borg "${script_dir}/consensus_sequence.py" "${script_dir}" ::: UBL*BORG*R2*.fastq.gz

# SERLOIN barcodes
parallel bash "${script_dir}/get_U6_barcodes.sh" {} serloin "${script_dir}/consensus_sequence.py" "${script_dir}" ::: UBL*SERLOIN*R2*.fastq.gz

# Random barcodes
parallel bash "${script_dir}/get_U6_barcodes.sh" {} random "${script_dir}/consensus_sequence.py" "${script_dir}" ::: UBL*50bp*R2*.fastq.gz

# Move results all into 1 parent directory
mkdir -p "${input_dir}/results"
rm -r 'UBL*BORG*R2*.fastq.gz_results' 'UBL*SERLOIN*R2*.fastq.gz_results' 'UBL*50bp*R2*.fastq.gz_results'
mv UBL*results results

# Make csv file of consensus sequences and calculate %GC content
bash consensus_gc_content.sh
mv barcode_consensus_sequences.csv results
mv trimmed_fastqs results

# Make multiqc report for mapping stats
cd results
multiqc .
