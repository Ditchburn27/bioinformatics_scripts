#! /usr/bin/bash

# bash script for TPM calculation using Stringtie

SAMPLES="RL3314 RL3315 RL3316 RL3317 RL3318 RL3319 RL3320 RL3321 RL3322 RL3323 RL3324 RL3325 RL3326 RL3327 RL3328 RL3329 RL3330 RL3331 RL3332 RL3333 RL3334 RL3335 RL3336 RL3337"

for SAMPLE in $SAMPLES; do
    stringtie -p 10 -G /home/lditchburn/working_data_02/rundata/reference/gencode.v19.annotation.gtf -o ${SAMPLE}-annot.gtf ${SAMPLE}.sorted.bam
done 