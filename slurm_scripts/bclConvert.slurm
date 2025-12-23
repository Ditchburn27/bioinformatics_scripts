#!/bin/bash --login
#SBATCH --job-name=Demux_BCL_convert
#SBATCH --partition=ll
#SBATCH --mem-per-cpu=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00   
#SBATCH --export=NONE
#SBATCH --mail-user=leighton.ditchburn@research.uwa.edu.au
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/group/ll005/lditchburn/SLURM_Logs/slurm_log-%x-%j.out


# To compile with the GNU toolchain
module purge
module load gcc/12.4.0 
module load Anaconda3/2024.06
conda activate /group/ll005/envs/cutntag


# Increase open file limit
ulimit -n 65535

# Check it
ulimit -n



# leave in, it lists the environment loaded by the modules
module list

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
if [ -z "$1" ]
then
      echo -ne "
      ################\n
      First argument is empty. Please provide path to Run folder.\n
      ################\n"
      exit 1
else
      RUNFOLDER=$1

      echo -ne "\nPath to Run folder is\n\t${RUNFOLDER}\n"
fi

if [ -z "$2" ]
then
      BARCODEMISMATCH=1	
      echo -ne "
      ################\n
      No barcode mismatch given for the second argument. Default is barcode mismatch = 1 \n
      ################\n"
else
      BARCODEMISMATCH=$2

      echo -ne "\nBarcode mismatch set to: ${BARCODEMISMATCH}\n"
fi

if [ -z "$3" ]
then
      SampleSheet="SampleSheet.csv"	
      echo -ne "
      ################\n
      No SampleSheet path given. Default is RunFolderPath/SampleSheet.csv \n
      ################\n"
else
      SampleSheet=$3

      echo -ne "\nSampleSheet is set to: ${SampleSheet}\n"
fi

if [ -z "$4" ]
then
      NoLaneSplit="true"	
      echo -ne "
      ################\n
      No Lane splitting option given. Defualt is true. \n
      ################\n"
else
      NoLaneSplit=$4

      echo -ne "\nno-lane-splitting is set to: ${NoLaneSplit}\n"
fi

JOBNAME=fastq_conversion
SCRATCH="/scratch/ll005/lditchburn/${JOBNAME}_${SLURM_JOBID}"
RESULTS="${RUNFOLDER}/fastqs_$(date +%Y-%m-%d_%H-%M-%S)"

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

# Create a Log folder for SLURM logs
mkdir ${RESULTS}/logs

###############################################
# Start of job
echo $JOBNAME job started at  `date`

################################################
# declare the name of the output file or log file
OUTPUT=${JOBNAME}_${SLURM_JOBID}.log

#############################################
# Run Demux

cd $SCRATCH/

bcl-convert --bcl-input-directory ${RUNFOLDER} \
--bcl-num-conversion-threads 12 \
--bcl-num-parallel-tiles 3 \
--bcl-num-compression-threads 4 \
--output-directory $SCRATCH/fastqs \
--bcl-sampleproject-subdirectories true \
--force \
--no-lane-splitting ${NoLaneSplit} \
--sample-sheet $SampleSheet

#############################################
# fastqc
module load fastqc
mkdir -p $SCRATCH/fastqc
fastqc -t $SLURM_CPUS_PER_TASK -o $SCRATCH/fastqc $SCRATCH/fastqs/*/*_R[12]_*fastq.gz

# multi fastqc
mkdir -p $SCRATCH/multiqc
multiqc -o $SCRATCH/multiqc $SCRATCH 


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
mv $SCRATCH/fastqs/* $RESULTS
mv $SCRATCH/fastqc $RESULTS
mv $SCRATCH/multiqc $RESULTS

# Delete SRATCH
rm -r $SCRATCH

echo $JOBNAME job finished at  `date`

# Copy the log file to Results
cp /group/ll005/lditchburn/SLURM_Logs/slurm_log-${SLURM_JOB_NAME}-${SLURM_JOBID}.out ${RESULTS}/logs

# Open up permissions
chown -R lditchburn:llusers ${RESULTS} 
chmod -R ug+rx ${RESULTS} 
