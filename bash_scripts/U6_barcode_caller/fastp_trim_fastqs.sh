#!/bin/bash

input_dir="$1"
sampleID="$2"
mkdir -p "${input_dir}/trimmed_fastqs"
output_folder="${input_dir}/trimmed_fastqs"
echo "Start Trimming..."
echo "Output folder: $output_folder"

trim() {
    local sample="$1"
    local output_folder="$2"  # Pass output_folder explicitly
    echo "Trimming sample: $sample"
    echo "Output folder inside function: $output_folder"
    fastp -i "$input_dir/${sample}"_R1*.gz -I "$input_dir/${sample}"_R2*.gz \
    --detect_adapter_for_pe -o "$output_folder/${sample}_R1.fastq.gz" \
    -O "$output_folder/${sample}_R2.fastq.gz" \
    -h "$output_folder/${sample}_fastp_report.html" \
    -j "$output_folder/${sample}_fastp_report.json" \
    -R "$output_folder/${sample}_fastp_report"

}

export -f trim

# Read sample IDs line by line and trim each sample in parallel
while IFS= read -r sample; do
    trim "$sample" "$output_folder"
done < "$sampleID"

echo "Done trimming..."