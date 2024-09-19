#!/bin/bash --login

#SBATCH --job-name=cellranger_RNA_count
#SBATCH --partition=peb
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

# Default chemistry value
CHEMISTRY="auto"

# Parse command line options
while getopts ":c:" opt; do
  case ${opt} in
    c )
      CHEMISTRY=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

# Take positional arguments (input_files, TRANSCRIPTOME, SAMPLE_FILE)
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
if [ -z "$3" ]; then
  echo -ne "
  ################\n
  Third argument is empty. Please provide sample file (sample.txt).\n
  ################\n"
  exit 1
else
  SAMPLE_FILE=$3
fi

# Create a unique directory in the SCRATCH directory for this job to run in.
if [ ! -d $SCRATCH ]; then 
  mkdir -p $SCRATCH 
fi 
echo "SCRATCH is $SCRATCH"

cp -r $input_files ${SCRATCH}/
cp $SAMPLE_FILE ${SCRATCH}

cd $SCRATCH

# Define FASTQ directory
FASTQ_DIR='fastqs'

# Read sample names from the sample.txt file and run cellranger count for each sample
while IFS= read -r SAMPLE_NAME; do
  echo "Running cellranger count for sample: $SAMPLE_NAME"
  
  cellranger count \
    --id="$SAMPLE_NAME" \
    --fastqs=$FASTQ_DIR \
    --sample="$SAMPLE_NAME" \
    --transcriptome=$TRANSCRIPTOME \
    --chemistry=$CHEMISTRY
done < "$SAMPLE_FILE"
