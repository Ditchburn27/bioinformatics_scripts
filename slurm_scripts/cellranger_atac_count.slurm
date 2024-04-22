#!/bin/bash --login

#SBATCH --job-name=cellranger_atac_count
#SBATCH --partition=peb
#SBATCH --mem-per-cpu=25G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --time=24:00:00   
#SBATCH --export=NONE
##SBATCH --mail-user=22720224@student.uwa.edu.au
##SBATCH --mail-type=BEGIN,END

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
JOBNAME=${SLURM_JOB_NAME}
SCRATCH=$MYSCRATCH/$JOBNAME/$SLURM_JOBID
REFERENCE=/group/ll005/reference/GRCh38_2020_atac

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
# create a sample ID file
ls fastqs/*/ > input_files.txt
sed 's/.*\(RL[0-9]\{4\}\).*/\1/' input_files.txt > sample_ids.txt
ls -d $PWD/fastqs/*/* > fastq_files.txt
fastq_file="fastq_files.txt"
sample_ids="sample_ids.txt"

###############################################
# Command to run cellranger ATAC count pipeline
command="cellranger-atac count --reference=$REFERENCE"

# Read each line from input.txt and output.txt simultaneously
while IFS= read -r fastq_dir && IFS= read -r sample_id <&3; do
    # Run the command for each pair of directories and sample_ids
    $command --id=$sample_id --fastqs=$fastq_dir
done < "$fastq_file" 3< "$sample_ids"
