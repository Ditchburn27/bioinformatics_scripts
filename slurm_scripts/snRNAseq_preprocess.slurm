#!/bin/bash --login
#SBATCH --job-name=snRNA-seq_preprocess
#SBATCH --partition=peb
#SBATCH --mem-per-cpu=10G
##SBATCH --reservation=training
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=4:00:00   
#SBATCH --export=NONE
#SBATCH --mail-user=22720224@student.uwa.edu.au
#SBATCH --mail-type=BEGIN,END

# To compile with the GNU toolchain
#module load gcc/9.4.0
module load Anaconda3/2020.11
conda activate /group/ll005/envs/scenicplus

# leave in, it lists the environment loaded by the modules
module list

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
if [ -z "$1" ]
then
      echo -ne "
      ################\n
      First argument is empty. Please provide path to 10X h5 matrix file.\n
      ################\n"
      exit 1
else
      FILE=$1


      echo -ne "\nPath to  path to 10X h5 matrix file is\n\t${FILE}\n"
fi

if [ -z "$2"]
then
     echo -ne "
     #################\n
     Second argument is empty. Please provide a name to prefix the output files.\n
     #################\n"
     exit 1
else
     OUTPUT_FILENAME=$2
     
     echo -ne "\nOutput files will be prefixed with \n\t${OUTPUT_FILENAME}\n"
fi


BASENAME=`basename $1`
JOBNAME=${OUTPUT_FILENAME}
SCRATCH=$MYSCRATCH/$JOBNAME/$SLURM_JOBID
RESULTS=$MYGROUP/$JOBNAME/$SLURM_JOBID
TMP=$MYSCRATCH/$JOBNAME/$SLURM_JOBID/TMP

###############################################
# Creates a unique directory in the SCRATCH directory for this job to run in.
if [ ! -d $SCRATCH ]; then 
    mkdir -p $SCRATCH
fi 
echo SCRATCH is $SCRATCH

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

#############################################
#   Copy input files to $SCRATCH
#   then change directory to $SCRATCH


cp $FILE $SCRATCH

cd $SCRATCH/


python3 /group/ll005/lditchburn/bioinformatics_scripts/snRNA-seq/snRNAseq_preprocess_v1.py $FILE -o${OUTPUT_FILENAME}

#############################################

# Get an idea of the directory structure and files
ls -lR * >> ${OUTPUT}

#############################################
#    $OUTPUT file to the unique results dir
# note this can be a copy or move  
mv  $OUTPUT ${RESULTS}

cd $HOME

###########################
# Clean up $SCRATCH 

# Save files that were created
mv $SCRATCH/* $RESULTS

# Delete SRATCH
rm -r $SCRATCH

echo $JOBNAME job finished at  `date`