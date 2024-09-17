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

JOBNAME=${SLURM_JOB_NAME}
SCRATCH=$MYSCRATCH/$JOBNAME/$SLURM_JOBID
################################################ 
# Take command line inputs
if [ -z "$1" ]; then
  echo -ne "
  ################\n
  First argument is empty. Please provide input file directory.\n
  ################\n"
  exit 1
else
  input_files=$1
fi
if [ -z "$2" ]; then
  echo -ne "
  ################\n
  Second argument is empty. Please provide reference file directory.\n
  ################\n"
  exit 1
else
  TRANSCRIPTOME=$2
fi

###############################################
# Create a unique directory in the SCRATCH directory for this job to run in.
if [ ! -d $SCRATCH ]; then 
  mkdir -p $SCRATCH 
fi 
echo SCRATCH is $SCRATCH

cp -r $input_files ${SCRATCH}/

cd $SCRATCH

###############################################
# Run cellranger count
FASTQ_DIR='fastqs'

# Extract the correct sample name pattern from the FASTQ file name
for file in $FASTQ_DIR/*.fastq.gz; do
  # Extract the sample name without path and file extension
  SAMPLE_NAME=$(basename "$file" | sed -E 's/(.*)_S[0-9]+_L[0-9]+_R[12]_001\.fastq\.gz/\1/' | sort | uniq)
  
  # Create a valid sample ID by replacing invalid characters with underscores
  SAMPLE_ID=$(echo "$SAMPLE_NAME" | sed 's/[^a-zA-Z0-9_-]/_/g')
  
  # Run cellranger count for each sample
  echo "Running cellranger count for sample: $SAMPLE_ID"
  
  cellranger count \
    --id=$SAMPLE_ID \
    --fastqs=$FASTQ_DIR \
    --sample=$SAMPLE_NAME \
    --transcriptome=$TRANSCRIPTOME
done
