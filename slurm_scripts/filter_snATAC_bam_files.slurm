#!/bin/bash --login

#SBATCH --job-name=bam_filter
#SBATCH --partition=peb
#SBATCH --mem-per-cpu=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --time=30:00:00   
#SBATCH --export=NONE
##SBATCH --mail-user=22720224@student.uwa.edu.au
##SBATCH --mail-type=BEGIN,END

# Start of job
echo $SLURM_JOB_NAME job started at  `date`

# To compile with the GNU toolchain
module load Anaconda3/2020.11
conda activate /group/ll005/envs/cutntag

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
JOBNAME=${SLURM_JOB_NAME}

# Define the directory containing BAM files
input_dir=$1

# Define the path to the blacklist BED file
blacklist_file="$SHARED/reference/bowtie2_hg38/blacklist/hg38-blacklist.v2.bed"

# Move to input_dir
cd $input_dir

# Run filtering
ls | parallel 'samtools view -b -f 0x2 -q 30 -h {}/outs/possorted_bam.bam > {}/outs/quality_filtered.bam'
ls | paralell 'samtools view -b -h -L $blacklist_file {}/outs/quality_filtered.bam -U {}/outs/blacklist_quality_filtered.bam'
ls | parallel 'samtools index {}/outs/fblacklist_quality_filtered.bam'

# End of job
echo $SLURM_JOB_NAME job finished at  `date`
