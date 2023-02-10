#!/bin/bash

# Script to merge sorted bam files and then calculate counts for genes

SECONDS=0

SAMPLES="RL3314 RL3315 RL3316 RL3317 RL3318 RL3319 RL3320 RL3321 RL3322 RL3323 RL3324 RL3325 RL3326 RL3327 RL3328 RL3329 RL3330 RL3331 RL3332 RL3333 RL3334 RL3335 RL3336 RL3337"

# STEP 1: merge bam files 
for SAMPLE in $SAMPLES; do
samtools merge merged_${SAMPLE}.bam ${SAMPLE}.bam ${SAMPLE}_v2.bam
done
echo "Bam files finished merging!"

# STEP 3: Run featureCounts - Quantification
# Gencode gtf for hg19 (2013-12-05) combined with gtf for ERCC spike-ins and transgenes
# -s 2 for TrueSeq mRNA strandard library prep
featureCounts -T 10 -s 2 -p -a /home/lditchburn/working_data_02/rundata/reference/hg19_hChr2_110_103_ERCC/gencodeV19_110_103_hChr2_ERCC.gtf -o extra_reads_counts_v5.txt merged_RL*.bam
echo "featureCounts finished running!"


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."