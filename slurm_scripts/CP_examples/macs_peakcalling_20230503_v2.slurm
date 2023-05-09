#!/bin/bash --login
#SBATCH --job-name=macs2
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
conda activate /group/peb002/conda_environments/macs2

# leave in, it lists the environment loaded by the modules
module list

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables

JOBNAME=macs2_peakcalling
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

cd $SCRATCH/

# Call ChIPseq peaks with Macs2

macs2 callpeak \
	-t /group/peb002/cpflueger/RL3810_P5C10_input/86841/RL3810_P5C10_input.bam \
	-c /group/peb002/cpflueger/RL3809_dCas9SunTag_HA_ChIP_P5C10/86840/RL3809_dCas9SunTag_HA_ChIP_P5C10.bam \
	 -g hs -f BAMPE -n 15xgRNA_P5C2_dCas9SunTag_ChIPseq \
	 -p 1e-06 --bdg --tempdir ~/scratch/tmp/ > 15xgRNA_P5C2_dCas9SunTag_ChIPseq_macs2.log 2>&1 &

wait

pigz *bdg

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

