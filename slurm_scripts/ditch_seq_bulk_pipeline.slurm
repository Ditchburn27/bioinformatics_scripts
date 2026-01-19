#!/bin/bash --login

#SBATCH --job-name=DITCH-seq_bulk
#SBATCH --partition=ll
#SBATCH --mem-per-cpu=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=/group/ll005/lditchburn/SLURM_Logs/slurm_log-%x-%j.out

################################################################################
# DITCH-seq Integrated Bulk Pipeline
# Author: Optimized for bulk reaction testing
# Features: UMI attachment, optional demultiplexing, full CUT&Tag analysis
################################################################################

echo "============================================================"
echo "DITCH-seq bulk pipeline started at $(date)"
echo "============================================================"

# Load modules
module purge
module load gcc
module load Anaconda3/2024.06
conda activate /group/ll005/envs/cutntag

################################################################################
# HARDCODED PATHS - EDIT THESE FOR YOUR SETUP
################################################################################

# Python scripts (edit these paths)
UMI_ATTACH_SCRIPT="/group/ll005/lditchburn/github_repos/DITCH-seq/bulk_script_folder/UMI_attach_EasySci_ATAC3.py"
DEMUX_SCRIPT="/group/ll005/lditchburn/github_repos/DITCH-seq/bulk_script_folder/demux_tn5.py"

# Default barcode pickle files
DEFAULT_RT_BARCODE_N5="/group/ll005/lditchburn/github_repos/DITCH-seq/script_folder/barcode_files/BC_N5.pickle2"
DEFAULT_RT_BARCODE_N7="/group/ll005/lditchburn/github_repos/DITCH-seq/script_folder/barcode_files/BC_N7.pickle2"
DEFAULT_LIGATION_BARCODE="/group/ll005/lditchburn/github_repos/DITCH-seq/script_folder/barcode_files/EasySci-ATAC_P5_BC.pickle2"

# Reference genome parameters
BOWTIE2_INDEX="/group/ll005/reference/bowtie2_mm10/mm10/mm10"
BLACKLIST="/group/ll005/reference/bowtie2_mm10/blacklist/mm10-blacklist.v2.bed"
GENOME_SIZE="mm"  # Change to "hs" for human

################################################################################
# INPUT PARAMETERS
################################################################################

# Required arguments
RAW_FASTQ_DIR=$1        # Directory containing raw FASTQ files (R1, R2, R3)
SAMPLE_LIST=$2          # File containing sample IDs (one per line, e.g., RL7238)
OUTPUT_BASE=$3          # Base output directory
MODE=$4                 # "singleplex" or "multiplex"

# Optional arguments (use defaults if not provided)
RT_BARCODE_N5=${5:-$DEFAULT_RT_BARCODE_N5}
RT_BARCODE_N7=${6:-$DEFAULT_RT_BARCODE_N7}
LIGATION_BARCODE=${7:-$DEFAULT_LIGATION_BARCODE}
TN5_CONFIG=${8:-""}     # Format: RL_NUMBER<tab>BARCODE_FILE_PATH

# Validate required inputs
if [[ -z "$RAW_FASTQ_DIR" ]] || [[ -z "$SAMPLE_LIST" ]] || [[ -z "$OUTPUT_BASE" ]] || [[ -z "$MODE" ]]; then
    echo "ERROR: Missing required arguments"
    echo ""
    echo "Usage: sbatch $0 <raw_fastq_dir> <sample_list> <output_dir> <mode> [RT_N5] [RT_N7] [LIGATION] [TN5_CONFIG]"
    echo ""
    echo "Required Arguments:"
    echo "  raw_fastq_dir  : Directory with raw FASTQ files"
    echo "  sample_list    : File with sample IDs (one per line)"
    echo "  output_dir     : Base output directory"
    echo "  mode          : 'singleplex' or 'multiplex'"
    echo ""
    echo "Simple Usage:"
    echo "  sbatch $0 raw_fastq/ samples.txt output/ singleplex"
    exit 1
fi

# Validate mode
if [[ "$MODE" != "singleplex" ]] && [[ "$MODE" != "multiplex" ]]; then
    echo "ERROR: Mode must be 'singleplex' or 'multiplex'"
    exit 1
fi

# Check if Python scripts exist
if [[ ! -f "$UMI_ATTACH_SCRIPT" ]]; then
    echo "ERROR: UMI attach script not found: $UMI_ATTACH_SCRIPT"
    echo "Please edit line 32 in this script with correct path"
    exit 1
fi

if [[ "$MODE" == "multiplex" ]] && [[ ! -f "$DEMUX_SCRIPT" ]]; then
    echo "ERROR: Demux script not found: $DEMUX_SCRIPT"
    echo "Please edit line 33 in this script with correct path"
    exit 1
fi

# Check if barcode files exist
if [[ ! -f "$RT_BARCODE_N5" ]]; then
    echo "ERROR: N5 barcode file not found: $RT_BARCODE_N5"
    exit 1
fi
if [[ ! -f "$RT_BARCODE_N7" ]]; then
    echo "ERROR: N7 barcode file not found: $RT_BARCODE_N7"
    exit 1
fi
if [[ ! -f "$LIGATION_BARCODE" ]]; then
    echo "ERROR: Ligation barcode file not found: $LIGATION_BARCODE"
    exit 1
fi

# Setup scratch directory
JOBNAME=${SLURM_JOB_NAME}
SCRATCH=$MYSCRATCH/$JOBNAME/$SLURM_JOBID

if [ ! -d $SCRATCH ]; then
    mkdir -p $SCRATCH
fi

echo "Working directory: $SCRATCH"
echo "Mode: $MODE"
echo "CPU threads: ${SLURM_CPUS_PER_TASK}"
echo "UMI script: $UMI_ATTACH_SCRIPT"
echo "N5 barcodes: $RT_BARCODE_N5"
echo "N7 barcodes: $RT_BARCODE_N7"
echo "Ligation barcodes: $LIGATION_BARCODE"

# Create directory structure
mkdir -p ${SCRATCH}/{raw_fastq,umi_attached,demuxed,fastqc,trimmed_fastq,sorted_bam,dedup_bam,fragment_bed,bigwig}
mkdir -p ${SCRATCH}/MACS2/{narrow,broad}
mkdir -p ${SCRATCH}/logs/{umi_attach,demux,fastp,bowtie2,dedup}
mkdir -p ${OUTPUT_BASE}

################################################################################
# STEP 0: Copy raw FASTQ files
################################################################################

echo ""
echo "============================================================"
echo "STEP 0: Copying raw FASTQ files"
echo "============================================================"

cp ${RAW_FASTQ_DIR}/*.fastq.gz ${SCRATCH}/raw_fastq/ 2>/dev/null || cp ${RAW_FASTQ_DIR}/*.fq.gz ${SCRATCH}/raw_fastq/ 2>/dev/null
cp $SAMPLE_LIST $SCRATCH

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to copy FASTQ files"
    exit 1
fi

cd $SCRATCH

################################################################################
# STEP 1: UMI Attachment
################################################################################

echo ""
echo "============================================================"
echo "STEP 1: UMI Attachment"
echo "============================================================"

python ${UMI_ATTACH_SCRIPT}     ${SCRATCH}/raw_fastq     ${SAMPLE_LIST}     ${SCRATCH}/umi_attached     ${RT_BARCODE_N5}     ${RT_BARCODE_N7}     ${LIGATION_BARCODE}     ${SLURM_CPUS_PER_TASK}     2>&1 | tee ${SCRATCH}/logs/umi_attach/umi_attach.log

if [ $? -ne 0 ]; then
    echo "ERROR: UMI attachment failed"
    exit 1
fi

FASTQ_SOURCE="${SCRATCH}/umi_attached"

################################################################################
# STEP 2: Demultiplexing by Tn5 index (multiplex mode only)
################################################################################

if [[ "$MODE" == "multiplex" ]]; then
    echo ""
    echo "============================================================"
    echo "STEP 2: Demultiplexing by Tn5 index"
    echo "============================================================"

    if [[ -n "$TN5_CONFIG" ]] && [[ -f "$TN5_CONFIG" ]]; then
        python ${DEMUX_SCRIPT}             ${FASTQ_SOURCE}             ${SAMPLE_LIST}             ${SCRATCH}/demuxed             ${TN5_CONFIG}             ${SLURM_CPUS_PER_TASK}             2>&1 | tee ${SCRATCH}/logs/demux/demux.log

        if [ $? -eq 0 ]; then
            FASTQ_SOURCE="${SCRATCH}/demuxed"
        else
            echo "WARNING: Demultiplexing failed, continuing with non-demuxed files"
        fi
    else
        echo "Warning: No Tn5 config file provided for multiplex mode"
    fi
else
    echo "Skipping demultiplexing (singleplex mode)"
fi

################################################################################
# STEP 3: FastQC on processed FASTQs
################################################################################

echo ""
echo "============================================================"
echo "STEP 3: FastQC"
echo "============================================================"

fastqc ${FASTQ_SOURCE}/*fastq.gz -o ${SCRATCH}/fastqc --threads ${SLURM_CPUS_PER_TASK}

################################################################################
# STEP 4: Adapter Trimming with fastp
################################################################################

echo ""
echo "============================================================"
echo "STEP 4: Adapter Trimming"
echo "============================================================"

ls ${FASTQ_SOURCE}/*R1*fastq.gz > ${SCRATCH}/read1_list.txt 2>/dev/null || ls ${FASTQ_SOURCE}/*R1*fq.gz > ${SCRATCH}/read1_list.txt
ls ${FASTQ_SOURCE}/*R2*fastq.gz > ${SCRATCH}/read2_list.txt 2>/dev/null || ls ${FASTQ_SOURCE}/*R2*fq.gz > ${SCRATCH}/read2_list.txt

export SLURM_CPUS_PER_TASK
parallel --jobs 4 --link     "fastp -w $((SLURM_CPUS_PER_TASK/4)) -x -y     -i {1} -I {2}     -o ${SCRATCH}/trimmed_fastq/{1//.fastq.gz}_trimmed.fastq.gz     -O ${SCRATCH}/trimmed_fastq/{2//.fastq.gz}_trimmed.fastq.gz     -h ${SCRATCH}/trimmed_fastq/{1//.fastq.gz}_fastp.html     -R ${SCRATCH}/trimmed_fastq/{1//.fastq.gz}_fastp_report     &> ${SCRATCH}/logs/fastp/{1//.fastq.gz}_fastp.log"     :::: ${SCRATCH}/read1_list.txt     :::: ${SCRATCH}/read2_list.txt

################################################################################
# STEP 5: Alignment with Bowtie2
################################################################################

echo ""
echo "============================================================"
echo "STEP 5: Bowtie2 Alignment"
echo "============================================================"

ls ${SCRATCH}/trimmed_fastq/*R1*trimmed.fastq.gz > ${SCRATCH}/trimmed_read1_list.txt
ls ${SCRATCH}/trimmed_fastq/*R2*trimmed.fastq.gz > ${SCRATCH}/trimmed_read2_list.txt

parallel --jobs 3 --link     "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant     --phred33 -I 10 -X 700 -p $((SLURM_CPUS_PER_TASK/3))     -x ${BOWTIE2_INDEX} -1 {1} -2 {2}     2> ${SCRATCH}/logs/bowtie2/{1//.fastq.gz}_bowtie2.txt |     samtools view -bS -F 4 - |     samtools sort -@ 4 -o ${SCRATCH}/sorted_bam/{1//.R1_001_trimmed.fastq.gz}.sorted.bam -"     :::: ${SCRATCH}/trimmed_read1_list.txt     :::: ${SCRATCH}/trimmed_read2_list.txt

parallel -j8 "samtools index {}" ::: ${SCRATCH}/sorted_bam/*.bam

################################################################################
# STEP 6: Remove Duplicates
################################################################################

echo ""
echo "============================================================"
echo "STEP 6: Removing Duplicates"
echo "============================================================"

parallel -j10     "picard MarkDuplicates -I {} -O ${SCRATCH}/dedup_bam/{/.}.rmDup.bam     --REMOVE_DUPLICATES true -M ${SCRATCH}/logs/dedup/{/.}_picard.txt"     ::: ${SCRATCH}/sorted_bam/*.bam

parallel -j8 "samtools index {}" ::: ${SCRATCH}/dedup_bam/*.bam

################################################################################
# STEP 7: Generate Fragment Files
################################################################################

echo ""
echo "============================================================"
echo "STEP 7: Generating Fragment Files"
echo "============================================================"

parallel -j10 "samtools sort -@ 2 -n {} -o ${SCRATCH}/fragment_bed/{/.}.name_sorted.bam"     ::: ${SCRATCH}/dedup_bam/*.bam

parallel -j10     "bedtools bamtobed -bedpe -i {} |     awk 'BEGIN{OFS=\"	\"}{if(\$1 !~ /chrM/)print \$1,\$2+4,\$6-5}' |     bedtools intersect -v -a - -b ${BLACKLIST} > ${SCRATCH}/fragment_bed/{/.}.bed"     ::: ${SCRATCH}/fragment_bed/*.name_sorted.bam

rm ${SCRATCH}/fragment_bed/*.name_sorted.bam

################################################################################
# STEP 8: Peak Calling with MACS2
################################################################################

echo ""
echo "============================================================"
echo "STEP 8: Peak Calling"
echo "============================================================"

parallel -j5     "macs2 callpeak -t {} -f BAMPE --keep-dup all --nolambda --nomodel     --gsize ${GENOME_SIZE} -q 1e-5 -n {/.} --outdir ${SCRATCH}/MACS2/narrow"     ::: ${SCRATCH}/dedup_bam/*.bam

parallel -j5     "macs2 callpeak -t {} -f BAMPE --keep-dup all --broad --nolambda --nomodel     --gsize ${GENOME_SIZE} -q 1e-5 -n {/.} --outdir ${SCRATCH}/MACS2/broad"     ::: ${SCRATCH}/dedup_bam/*.bam

################################################################################
# STEP 9: Generate BigWig Files
################################################################################

echo ""
echo "============================================================"
echo "STEP 9: Generating BigWig Files"
echo "============================================================"

parallel -j4     "bamCoverage -p $((SLURM_CPUS_PER_TASK/4)) -b {}     --normalizeUsing CPM --blackListFileName ${BLACKLIST}     -o ${SCRATCH}/bigwig/{/.}.bw"     ::: ${SCRATCH}/dedup_bam/*.bam

################################################################################
# STEP 10: QC and Summary Statistics
################################################################################

echo ""
echo "============================================================"
echo "STEP 10: Quality Control and Statistics"
echo "============================================================"

multiqc ${SCRATCH} -o ${SCRATCH} --filename ditch_seq_multiqc_report

if command -v seqkit &> /dev/null; then
    seqkit stats -j ${SLURM_CPUS_PER_TASK} ${FASTQ_SOURCE}/*.fastq.gz > ${SCRATCH}/processed_fastq_stats.txt 2>/dev/null || true
    seqkit stats -j ${SLURM_CPUS_PER_TASK} ${SCRATCH}/trimmed_fastq/*.fastq.gz > ${SCRATCH}/trimmed_fastq_stats.txt 2>/dev/null || true
fi

echo -e "BAM_file	Unique_fragments" > ${SCRATCH}/bam_stats_summary.txt
for bam in ${SCRATCH}/dedup_bam/*.bam; do
    count=$(samtools view -@ 4 -F 0x40 $bam | cut -f1 | sort -u | wc -l)
    echo -e "$bam	$count" >> ${SCRATCH}/bam_stats_summary.txt
done

################################################################################
# STEP 11: Copy Results to Output Directory
################################################################################

echo ""
echo "============================================================"
echo "STEP 11: Copying Results"
echo "============================================================"

rsync -av --progress ${SCRATCH}/ ${OUTPUT_BASE}/

# Copy the log file to Results
cp /group/ll005/lditchburn/SLURM_Logs/slurm_log-${SLURM_JOB_NAME}-${SLURM_JOBID}.out ${OUTPUT_BASE}/

################################################################################
# COMPLETION
################################################################################

echo ""
echo "============================================================"
echo "DITCH-seq bulk pipeline completed at $(date)"
echo "Results saved to: ${OUTPUT_BASE}"
echo "============================================================"

rm -rf ${SCRATCH}