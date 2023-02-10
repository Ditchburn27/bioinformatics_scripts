#! /usr/bin/bash

# bash script to convert sam files to bam files using samtools 

SAMPLES="RL3314 RL3315 RL3316 RL3317 RL3318 RL3319 RL3320 RL3321 RL3322 RL3323 RL3324 RL3325 RL3326 RL3327 RL3328 RL3329 RL3330 RL3331 RL3332 RL3333 RL3334 RL3335 RL3336 RL3337"

for SAMPLE in $SAMPLES; do
    samtools view -bS ${SAMPLE}*.sam > ${SAMPLE}.bam
    samtools sort ${SAMPLE}.bam -o ${SAMPLE}.sorted.bam
    samtools index ${SAMPLE}.sorted.bam
done