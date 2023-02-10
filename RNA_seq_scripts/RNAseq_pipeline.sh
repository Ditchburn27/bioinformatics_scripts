#!/bin/bash

# RNA-seq qc, alignment and feature counts table pipeline
# Updated gtf file for featureCounts 28-09-22

SECONDS=0

SAMPLES="RL3314 RL3315 RL3316 RL3317 RL3318 RL3319 RL3320 RL3321 RL3322 RL3323 RL3324 RL3325 RL3326 RL3327 RL3328 RL3329 RL3330 RL3331 RL3332 RL3333 RL3334 RL3335 RL3336 RL3337"

# STEP 1: Run fastqc
for SAMPLE in $SAMPLES; do
    fastqc ${SAMPLE}_2022_10_02*R1_001.fastq.gz -o /home/lditchburn/working_data_02/Neuronal_Maturation/
    fastqc ${SAMPLE}_2022_10_02*R2_001.fastq.gz -o /home/lditchburn/working_data_02/Neuronal_Maturation/
done
echo "fastqc finished running!"

# STEP 2: Run HISAT2 alignment
# Pipe output to samtools to create a sorted bam file
for SAMPLE in $SAMPLES; do 
    hisat2 -p10 --rna-strandness 'RF' --dta -x /home/lditchburn/working_data_02/rundata/reference/hg19_hChr2_110_103_ERCC/hg19_hChr2_110_103_ERCC_hisat2 -1 ${SAMPLE}_2022_10_02*R1_001.fastq.gz -2 ${SAMPLE}_2022_10_02*R2_001.fastq.gz | samtools sort -o rf_aligned/${SAMPLE}_v2.bam
done 
echo "HISAT2 finished running!"

# Change directory to location of bam files
cd rf_aligned

# STEP 3: Run featureCounts - Quantification
# Gencode gtf for hg19 (2013-12-05) combined with gtf for ERCC spike-ins and transgenes
# -s 2 for TrueSeq mRNA strandard library prep
featureCounts -T 10 -s 2 -p -a /home/lditchburn/working_data_02/rundata/reference/hg19_hChr2_110_103_ERCC/gencodeV19_110_103_hChr2_ERCC.gtf -o extra_reads_counts_v4.txt RL*.bam
echo "featureCounts finished running!"


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."