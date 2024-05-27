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

################################################
## Filter bam files from cellranger output
# Define the directory containing BAM files
input_dir=$1

# Define the path to the blacklist BED file
blacklist_file="$SHARED/reference/bowtie2_hg38/blacklist/hg38-blacklist.v2.bed"

# Move to input_dir
cd $input_dir

# Run filtering
echo "Starting bam filtering"
ls | parallel 'samtools view -b -f 0x2 -q 30 -h {}/outs/possorted_bam.bam > {}/outs/quality_filtered.bam'

# Index quality filter bam files
ls | parallel 'samtools index {}/outs/quality_filtered.bam'

echo "Finished filtering bam files"

################################################
## Create Fragment Files for downstream analysis
# sinto to convert bam into fragments file
echo "Starting to create fragment files"
ls | parallel 'sinto fragments -b {}/outs/quality_filtered.bam \
-f {}/outs/filtered_fragments.bed \
--barcode_regex "[^:]*" \
--collapse_within'

# sort fragment files
ls | parallel 'sort -k 1,1 -k2,2n {}/outs/filtered_fragments.bed > {}/outs/filtered_fragments.sorted.bed'

# compress sorted fragment files
ls | parallel 'gzip {}/outs/filtered_fragments.sorted.bed'
echo "Finished making fragment files"
################################################
## Create h5ad anndata files from fragment files
# with snapATAC2 for downstream analysis
conda activate /group/ll005/envs/snapatac
echo "Start creating anndata object h5ad files"

python /group/ll005/lditchburn/bioinformatics_scripts/snATAC_analysis_scripts/ 1_snapATAC_make_h5ad_files.py $input_dir

echo "Finished creating anndata object h5ad files"
# End of job
echo $SLURM_JOB_NAME job finished at  `date`
