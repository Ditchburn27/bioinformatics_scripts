#!/bin/bash --login
#SBATCH --job-name=scvi_model_training
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --export=NONE
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --mail-user=22720224@student.uwa.edu.au
#SBATCH --mail-type=BEGIN,END

# To compile with the GNU toolchain
#module load gcc/9.4.0
module load Anaconda3/2020.11
module load cuda/11.6 
conda activate /group/ll005/envs/scarches

# leave in, it lists the environment loaded by the modules
module list

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables

JOBNAME=scvi_model_training
RESULTS=$MYGROUP/$JOBNAME/$SLURM_JOBID
TMP=$MYSCRATCH/$JOBNAME/$SLURM_JOBID/TMP

###############################################
# Creates a unique directory in your GROUP directory for the results of this job
if [ ! -d $RESULTS ]; then 
    mkdir -p $RESULTS 
fi
echo the results directory is $RESULTS
###############################################
# Creates a unique TMP directory on SCRATCH

if [ ! -d $TMP ]; then
    mkdir -p $TMP
    export TMPDIR=$TMP
fi
echo the temp directory is $TMP and TMPDIR is $TMPDIR


###############################################
# Start of job
echo $JOBNAME job started at  `date`

################################################
# declare the name of the output file or log file
OUTPUT=${JOBNAME}_${SLURM_JOBID}.log


python3 scvi_model_training.py



#############################################

# Get an idea of the directory structure and files
ls -lR * >> ${OUTPUT}


#############################################
#    $OUTPUT file to the unique results dir
# note this can be a copy or move  
mv  $OUTPUT ${RESULTS}

cd $HOME

############################################
echo $JOBNAME job finished at  `date`

