#!/bin/bash --login

#SBATCH --job-name=model_training
#SBATCH --partition=gpu
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:1
#SBATCH --time=8:00:00   
#SBATCH --export=NONE
##SBATCH --mail-user=22720224@student.uwa.edu.au
##SBATCH --mail-type=BEGIN,END

set -eu -o pipefail -o verbose

# Start of job
echo $SLURM_JOB_NAME job started at  `date`

# To compile with the GNU toolchain
module load Anaconda3/2020.11
conda activate /group/ll005/envs/rapid_scrna

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
JOBNAME=${SLURM_JOB_NAME}

if [ -z "$1" ]
then
      echo -ne "
      ################\n
      First argument is empty. Please provide path to python script.\n
      ################\n"
      exit 1
else
      PYTHON_SCRIPT=$1


      echo -ne "\nPath to python script is\n\t${PYTHON_SCRIPT}\n"
fi

python $PYTHON_SCRIPT