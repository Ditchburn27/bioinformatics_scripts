#!/bin/bash --login

#SBATCH --job-name=cellranger_RNA_count
#SBATCH --partition=work
#SBATCH --mem-per-cpu=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00   
#SBATCH --export=NONE
##SBATCH --mail-user=22720224@student.uwa.edu.au
##SBATCH --mail-type=BEGIN,END

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
JOBNAME=${SLURM_JOB_NAME}
SCRATCH=$MYSCRATCH/$JOBNAME/$SLURM_JOBID
################################################ 
#Take command line inputs
if [ -z "$1" ]
then
      echo -ne "
      ################\n
      First argument is empty. Please provide input file directory.\n
      ################\n"
      exit 1
else
      input_files=$1
fi
if [ -z "$2" ]
then
      echo -ne "
      ################\n
      Second argument is empty. Please provide reference file directory.\n
      ################\n"
      exit 1
else
      TRANSCRIPTOME=$2
fi

###############################################
# Creates a unique directory in the SCRATCH directory for this job to run in.
if [ ! -d $SCRATCH ]; then 
    mkdir -p $SCRATCH 
fi 
echo SCRATCH is $SCRATCH

mkdir -p ${SCRATCH}/fastqs
cp -r $input_files ${SCRATCH}/fastqs/

cd $SCRATCH

###############################################
# Run cellranger count
FASTQ_DIR='fastqs'
# Loop through each unique sample name in the FASTQ file names
for SAMPLE_NAME in $(ls $FASTQ_DIR/*.fastq.gz | sed -E 's/.*\/([^_]+).*/\1/' | sort | uniq); do
  # Run cellranger count for each sample
  cellranger count \
    --id=$SAMPLE_NAME \
    --fastqs=$FASTQ_DIR \
    --sample=$SAMPLE_NAME \
    --transcriptome=$TRANSCRIPTOME
done
