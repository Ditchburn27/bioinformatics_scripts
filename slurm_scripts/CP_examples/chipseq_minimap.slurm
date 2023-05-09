#!/bin/bash --login
#SBATCH --job-name=chipseq
#SBATCH --partition=peb
#SBATCH --mem-per-cpu=10G
##SBATCH --reservation=training
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00   
#SBATCH --export=NONE
#SBATCH --mail-user=christian.pflueger@uwa.edu.au
#SBATCH --mail-type=BEGIN,END

# To compile with the GNU toolchain
#module load gcc/9.4.0
module load Anaconda3/2020.11
conda activate /group/peb002/conda_environments/bioinfo

# leave in, it lists the environment loaded by the modules
module list

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
if [ -z "$1" ]
then
      echo -ne "
      ################\n
      First argument is empty. Please provide path to read 1 fastq.\n
      ################\n"
      exit 1
else
      R1=$1


      echo -ne "\nPath to  path to read 1 fastq is\n\t${R1}\n"
fi

BASENAME=`basename $R1`
JOBNAME=${BASENAME%%_S[0-9]*_*}
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


cp $1 $2 $SCRATCH

cd $SCRATCH/


bash /home/cpflueger/scripts/chip-seq/ChIPseq_fastp_minimap2pipeline_PE_kaya.sh \
	$1 \
	$2 \
	$3 \
	$4 \
	$5 \
	$TMP

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

# Delete TMP folder
rm -r $TMP

# Save files that were created
mv $SCRATCH/* $RESULTS

# Delete SRATCH
rm -r $SCRATCH

echo $JOBNAME job finished at  `date`

