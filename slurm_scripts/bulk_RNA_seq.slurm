#!/bin/bash --login

#SBATCH --job-name=Bulk_RNAseq
#SBATCH --partition=peb
#SBATCH --mem-per-cpu=30G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00   
#SBATCH --export=NONE
##SBATCH --mail-user=22720224@student.uwa.edu.au
##SBATCH --mail-type=BEGIN,END

# Start of job
echo $SLURM_JOB_NAME job started at  `date`

# To compile with the GNU toolchain
# To compile with the GNU toolchain
module load Anaconda3/2023.07
conda activate rnaseq

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
JOBNAME=${SLURM_JOB_NAME}
SCRATCH=$MYSCRATCH/$JOBNAME/$SLURM_JOBID


# take input directory with fastqs
FASTQ=$1
INDEX=$2
GTF_file=$3
###############################################
# Creates a unique directory in the SCRATCH directory for this job to run in.
if [ ! -d $SCRATCH ]; then 
    mkdir -p $SCRATCH 
fi 
echo SCRATCH is $SCRATCH

#############################################
#   Copy input files to $SCRATCH
#   then change directory to $SCRATCH

mkdir -p ${SCRATCH}/{fastq,fastqc,trimmed_fastq,sorted_bam,featureCounts}
mkdir -p ${SCRATCH}/logs/{fastp,hisat2}

### COPY files to scratch
cp ${FASTQ}/*gz ${SCRATCH}/fastq

## IMPORTANT: change directory to scratch

cd $SCRATCH
###################
#1. FastqQC

fastqc fastq/*fastq.gz -o ${SCRATCH}/fastqc --threads 20

#2.get the read1 and read2 name and put in a file
ls fastq/*_R1*fastq.gz >> read1
ls fastq/*_R2*fastq.gz >> read2

#3.Adapter Trimming with fastp

for i in $(cat read1)
do 
	read j
	echo $i $j
	fastp_R1=$(basename ${i/.fastq.gz}_trimmed.fastq.gz)
	fastp_R2=$(basename ${j/.fastq.gz}_trimmed.fastq.gz)
	html=$(basename ${i%_R1_001.fastq.gz}_fastp_report.html)
	report=$(basename ${i%_R1_001.fastq.gz}_fastp_report)
	fastp_log=$(basename ${i%_R1_001.fastq.gz}_fastp.out)
	fastp -w 20 -x -y \
	-i $i \
	-I $j \
	-o ${SCRATCH}/trimmed_fastq/${fastp_R1} \
	-O ${SCRATCH}/trimmed_fastq/${fastp_R2} \
	-h ${SCRATCH}/trimmed_fastq/${html} \
	-R ${SCRATCH}/trimmed_fastq/${report} &> ${SCRATCH}/logs/fastp/${fastp_log}
done < read2
#4.get the read1 and read2 name of trimmed fastq and put in a file

ls trimmed_fastq/*_R1*fastq.gz >> trimmed_read1
ls trimmed_fastq/*_R2*fastq.gz >> trimmed_read2

#5. map with HiSat2
for i in $(cat trimmed_read1)
do 
	read j
	echo $i $j
    hisat2_log=$(basename ${i%_R1_001_trimmed.fastq.gz}_hisat2.txt)
	sorted_bam=$(basename ${i%_R1_001_trimmed.fastq.gz}.sorted.bam)
    hisat2 -p 20 --rna-strandness 'RF' --dta -x $INDEX -1 $i -2 $j --summary-file logs/hisat2/${hisat2_log}| samtools view -bS -F 4 - | samtools sort -o sorted_bam/${sorted_bam} - 
    done < trimmed_read2
parallel -j5 "samtools index {} {.}.bam.bai" ::: sorted_bam/*bam

#6. featureCounts
featureCounts -T 20 -s 2 -p -a ${GTF_file} -o featureCounts/counts.txt sorted_bam/*.sorted.bam

#7. multiqc 
conda activate cutntag
multiqc .

echo $JOBNAME job finished at  `date`