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
REFERENCE=/group/ll005/reference/GRCh38-2020

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
# Run cellranger count
cellranger count --id=RL2290 --fastqs=fastqs/RL2290_2021_03_11_Organoid_12mo_DigitoninPostSort_10xTEST_snRNAseq_v2 --sample=RL2290_2021_03_11_Organoid_12mo_DigitoninPostSort_10xTEST_snRNAseq --transcriptome=$REFERENCE