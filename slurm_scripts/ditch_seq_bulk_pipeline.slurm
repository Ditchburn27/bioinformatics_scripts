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
# HARDCODED BARCODE PATHS - EDIT THESE FOR YOUR LAB
################################################################################

# Default barcode pickle files (edit these paths for your setup)
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
    echo "  raw_fastq_dir  : Directory with raw FASTQ files (R1.fastq.gz, R2.fastq.gz, R3.fastq.gz)"
    echo "  sample_list    : File with sample IDs (e.g., RL7238, one per line)"
    echo "  output_dir     : Base output directory"
    echo "  mode          : 'singleplex' or 'multiplex'"
    echo ""
    echo "Optional Arguments (will use defaults if not provided):"
    echo "  RT_N5         : N5 RT barcode pickle file"
    echo "                  Default: $DEFAULT_RT_BARCODE_N5"
    echo "  RT_N7         : N7 RT barcode pickle file"
    echo "                  Default: $DEFAULT_RT_BARCODE_N7"
    echo "  LIGATION      : Ligation barcode pickle file"
    echo "                  Default: $DEFAULT_LIGATION_BARCODE"
    echo "  TN5_CONFIG    : TSV file mapping RL numbers to Tn5 barcode files (multiplex only)"
    echo ""
    echo "Simple Usage Examples:"
    echo "  # Singleplex (using default barcodes)"
    echo "  sbatch $0 raw_fastq/ samples.txt output/ singleplex"
    echo ""
    echo "  # Multiplex (using default barcodes + Tn5 config)"
    echo "  sbatch $0 raw_fastq/ samples.txt output/ multiplex "" "" "" tn5_config.txt"
    echo ""
    echo "  # Override barcode files"
    echo "  sbatch $0 raw_fastq/ samples.txt output/ singleplex /path/to/N5.pkl /path/to/N7.pkl /path/to/lig.pkl"
    exit 1
fi

# Validate mode
if [[ "$MODE" != "singleplex" ]] && [[ "$MODE" != "multiplex" ]]; then
    echo "ERROR: Mode must be 'singleplex' or 'multiplex'"
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
cd $SCRATCH

################################################################################
# STEP 1: UMI Attachment
################################################################################

echo ""
echo "============================================================"
echo "STEP 1: UMI Attachment"
echo "============================================================"

# Create inline Python script for UMI attachment
cat > ${SCRATCH}/umi_attach.py << 'EOF'
import sys
import gzip
import pickle
from multiprocessing import Pool
from functools import partial

def attach_umi(sample, input_folder, output_folder, RT_N5, RT_N7, lig_bc):
    """Attach UMI barcodes to FASTQ headers with optimized I/O"""
    try:
        # Input files
        r1_in = f"{input_folder}/{sample}.R1.fastq.gz"
        r2_in = f"{input_folder}/{sample}.R3.fastq.gz" 
        r3_in = f"{input_folder}/{sample}.R2.fastq.gz"

        # Output files
        r1_out = f"{output_folder}/{sample}.R1.fastq.gz"
        r2_out = f"{output_folder}/{sample}.R2.fastq.gz"

        stats = {'total': 0, 'lig_match': 0, 'filtered': 0}

        with gzip.open(r1_in, 'rt') as f1,              gzip.open(r2_in, 'rt') as f2,              gzip.open(r3_in, 'rt') as f3,              gzip.open(r1_out, 'wt') as out1,              gzip.open(r2_out, 'wt') as out2:

            while True:
                # Read FASTQ records
                try:
                    h1, s1, p1, q1 = [next(f1).rstrip() for _ in range(4)]
                    h2, s2, p2, q2 = [next(f2).rstrip() for _ in range(4)]
                    h3, s3, p3, q3 = [next(f3).rstrip() for _ in range(4)]
                except StopIteration:
                    break

                stats['total'] += 1

                # Check ligation barcode (first 10bp of R2/I5)
                lig_bc_seq = s3[:10]
                if lig_bc_seq not in lig_bc:
                    continue

                stats['lig_match'] += 1
                lig_match = lig_bc[lig_bc_seq]

                # Check N5 barcode (first 6bp of R1)
                bc1_seq = s1[:6]
                if bc1_seq not in RT_N5:
                    continue
                bc1_match = RT_N5[bc1_seq]

                # Check N7 barcode (first 10bp of R3)
                bc2_seq = s2[:10]
                if bc2_seq not in RT_N7:
                    continue
                bc2_match = RT_N7[bc2_seq]

                stats['filtered'] += 1

                # Create new header with barcodes
                new_header = f"@{lig_match}{bc1_match}{bc2_match},{h1[1:]}"

                # Write R1 (trim first 25bp)
                out1.write(f"{new_header}
{s1[25:]}
{p1}
{q1[25:]}
")

                # Write R2 (trim first 29bp)  
                new_header2 = f"@{lig_match}{bc1_match}{bc2_match},{h2[1:]}"
                out2.write(f"{new_header2}
{s2[29:]}
{p2}
{q2[29:]}
")

        filter_rate = stats['filtered'] / stats['total'] * 100 if stats['total'] > 0 else 0
        print(f"{sample}: Total={stats['total']}, Lig_matched={stats['lig_match']}, "
              f"Filtered={stats['filtered']} ({filter_rate:.2f}%)")

        return stats

    except Exception as e:
        print(f"ERROR processing {sample}: {e}")
        return None

def main():
    input_dir, sample_list, output_dir, bc_n5, bc_n7, bc_lig, cores = sys.argv[1:]

    # Load barcodes
    with open(bc_n5, 'rb') as f:
        RT_N5 = pickle.load(f)
    with open(bc_n7, 'rb') as f:
        RT_N7 = pickle.load(f)
    with open(bc_lig, 'rb') as f:
        lig_bc = pickle.load(f)

    # Load samples
    with open(sample_list) as f:
        samples = [line.strip() for line in f if line.strip()]

    print(f"Processing {len(samples)} samples with {cores} cores...")

    # Parallel processing
    func = partial(attach_umi, input_folder=input_dir, output_folder=output_dir,
                   RT_N5=RT_N5, RT_N7=RT_N7, lig_bc=lig_bc)

    with Pool(int(cores)) as pool:
        results = pool.map(func, samples)

    print("
UMI attachment complete!")

if __name__ == "__main__":
    main()
EOF

python ${SCRATCH}/umi_attach.py     ${SCRATCH}/raw_fastq     ${SAMPLE_LIST}     ${SCRATCH}/umi_attached     ${RT_BARCODE_N5}     ${RT_BARCODE_N7}     ${LIGATION_BARCODE}     ${SLURM_CPUS_PER_TASK}     2>&1 | tee ${SCRATCH}/logs/umi_attach/umi_attach.log

FASTQ_SOURCE="${SCRATCH}/umi_attached"

################################################################################
# STEP 2: Demultiplexing by Tn5 index (multiplex mode only)
################################################################################

if [[ "$MODE" == "multiplex" ]]; then
    echo ""
    echo "============================================================"
    echo "STEP 2: Demultiplexing by Tn5 index"
    echo "============================================================"

    # Create inline Python script for demultiplexing
    cat > ${SCRATCH}/demux_tn5.py << 'EOF'
import sys
import gzip
import os
from collections import defaultdict
from multiprocessing import Pool
from functools import partial

def demux_sample(sample, input_dir, output_dir, barcode_dict, barcode_pos=10, barcode_len=6):
    """Demultiplex by Tn5 index in header"""

    r1_file = f"{input_dir}/{sample}.R1.fastq.gz"
    r2_file = f"{input_dir}/{sample}.R2.fastq.gz"

    if not os.path.exists(r1_file):
        print(f"Warning: {r1_file} not found")
        return None

    os.makedirs(output_dir, exist_ok=True)

    # File handles (created on-demand)
    handles_r1 = {}
    handles_r2 = {}
    unknown_r1 = unknown_r2 = None

    stats = defaultdict(int)
    total = 0

    with gzip.open(r1_file, 'rt') as f1, gzip.open(r2_file, 'rt') as f2:
        while True:
            try:
                h1, s1, p1, q1 = [next(f1).rstrip() for _ in range(4)]
                h2, s2, p2, q2 = [next(f2).rstrip() for _ in range(4)]
            except StopIteration:
                break

            total += 1

            # Extract Tn5 barcode from header
            if h1.startswith('@'):
                header_bc = h1[1:].split(',')[0]
                if len(header_bc) >= barcode_pos + barcode_len:
                    tn5_bc = header_bc[barcode_pos:barcode_pos + barcode_len]

                    if tn5_bc in barcode_dict:
                        stats[tn5_bc] += 1

                        # Create handles if needed
                        if tn5_bc not in handles_r1:
                            bc_name = barcode_dict[tn5_bc]
                            handles_r1[tn5_bc] = gzip.open(
                                f"{output_dir}/{sample}_{bc_name}_R1.fastq.gz", 'wt')
                            handles_r2[tn5_bc] = gzip.open(
                                f"{output_dir}/{sample}_{bc_name}_R2.fastq.gz", 'wt')

                        # Write reads
                        handles_r1[tn5_bc].write(f"{h1}
{s1}
{p1}
{q1}
")
                        handles_r2[tn5_bc].write(f"{h2}
{s2}
{p2}
{q2}
")
                        continue

            # Unknown barcode
            stats['unknown'] += 1
            if unknown_r1 is None:
                unknown_r1 = gzip.open(f"{output_dir}/{sample}_unknown_R1.fastq.gz", 'wt')
                unknown_r2 = gzip.open(f"{output_dir}/{sample}_unknown_R2.fastq.gz", 'wt')
            unknown_r1.write(f"{h1}
{s1}
{p1}
{q1}
")
            unknown_r2.write(f"{h2}
{s2}
{p2}
{q2}
")

    # Close handles
    for h in handles_r1.values():
        h.close()
    for h in handles_r2.values():
        h.close()
    if unknown_r1:
        unknown_r1.close()
        unknown_r2.close()

    # Print stats
    print(f"
{'='*60}")
    print(f"Sample: {sample}")
    print(f"Total reads: {total}")
    for bc in sorted(handles_r1.keys()):
        pct = stats[bc] / total * 100 if total > 0 else 0
        print(f"  {barcode_dict[bc]} ({bc}): {stats[bc]} ({pct:.2f}%)")
    if stats['unknown'] > 0:
        pct = stats['unknown'] / total * 100 if total > 0 else 0
        print(f"  Unknown: {stats['unknown']} ({pct:.2f}%)")
    print('='*60)

def main():
    input_dir, sample_list_file, output_dir, tn5_config, cores = sys.argv[1:]

    # Load samples
    with open(sample_list_file) as f:
        samples = [line.strip() for line in f if line.strip()]

    # Load Tn5 config (maps RL numbers to barcode files)
    rl_to_barcodes = {}
    if os.path.exists(tn5_config):
        with open(tn5_config) as f:
            for line in f:
                if line.strip():
                    rl, bc_file = line.strip().split('	')
                    # Load barcode file
                    bc_dict = {}
                    with open(bc_file) as bf:
                        for bcline in bf:
                            if bcline.strip():
                                parts = bcline.strip().split('	')
                                if len(parts) == 2:
                                    bc_dict[parts[0]] = parts[1]
                    rl_to_barcodes[rl] = bc_dict
                    print(f"Loaded {len(bc_dict)} barcodes for {rl}")

    # Default barcode dict if no config
    default_barcodes = None
    if not rl_to_barcodes:
        # Use a default if no config provided (will process all with same barcodes)
        print("Warning: No Tn5 config provided, skipping demux")
        return

    # Process each sample
    for sample in samples:
        # Find matching barcode dict
        barcode_dict = None
        for rl_prefix, bc_dict in rl_to_barcodes.items():
            if sample.startswith(rl_prefix):
                barcode_dict = bc_dict
                break

        if barcode_dict:
            demux_sample(sample, input_dir, output_dir, barcode_dict)
        else:
            print(f"Warning: No Tn5 barcodes found for {sample}")

if __name__ == "__main__":
    main()
EOF

    if [[ -n "$TN5_CONFIG" ]] && [[ -f "$TN5_CONFIG" ]]; then
        python ${SCRATCH}/demux_tn5.py             ${FASTQ_SOURCE}             ${SAMPLE_LIST}             ${SCRATCH}/demuxed             ${TN5_CONFIG}             ${SLURM_CPUS_PER_TASK}             2>&1 | tee ${SCRATCH}/logs/demux/demux.log

        FASTQ_SOURCE="${SCRATCH}/demuxed"
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

# Generate read pairs list
ls ${FASTQ_SOURCE}/*R1*fastq.gz > ${SCRATCH}/read1_list.txt
ls ${FASTQ_SOURCE}/*R2*fastq.gz > ${SCRATCH}/read2_list.txt

# Parallel trimming with fastp
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

# Parallel alignment (3 jobs to balance CPU usage)
parallel --jobs 3 --link     "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant     --phred33 -I 10 -X 700 -p $((SLURM_CPUS_PER_TASK/3))     -x ${BOWTIE2_INDEX} -1 {1} -2 {2}     2> ${SCRATCH}/logs/bowtie2/{1//.fastq.gz}_bowtie2.txt |     samtools view -bS -F 4 - |     samtools sort -@ 4 -o ${SCRATCH}/sorted_bam/{1//.R1_001_trimmed.fastq.gz}.sorted.bam -"     :::: ${SCRATCH}/trimmed_read1_list.txt     :::: ${SCRATCH}/trimmed_read2_list.txt

# Index BAM files
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

# Sort by name for bedpe conversion
parallel -j10 "samtools sort -@ 2 -n {} -o ${SCRATCH}/fragment_bed/{/.}.name_sorted.bam"     ::: ${SCRATCH}/dedup_bam/*.bam

# Convert to fragments (Tn5 shift: +4bp left, -5bp right)
parallel -j10     "bedtools bamtobed -bedpe -i {} |     awk 'BEGIN{OFS="	"}{if(\$1 !~ /chrM/)print \$1,\$2+4,\$6-5}' |     bedtools intersect -v -a - -b ${BLACKLIST} > ${SCRATCH}/fragment_bed/{/.}.bed"     ::: ${SCRATCH}/fragment_bed/*.name_sorted.bam

rm ${SCRATCH}/fragment_bed/*.name_sorted.bam

################################################################################
# STEP 8: Peak Calling with MACS2
################################################################################

echo ""
echo "============================================================"
echo "STEP 8: Peak Calling"
echo "============================================================"

# Narrow peaks
parallel -j5     "macs2 callpeak -t {} -f BAMPE --keep-dup all --nolambda --nomodel     --gsize ${GENOME_SIZE} -q 1e-5 -n {/.} --outdir ${SCRATCH}/MACS2/narrow"     ::: ${SCRATCH}/dedup_bam/*.bam

# Broad peaks
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

# MultiQC report
multiqc ${SCRATCH} -o ${SCRATCH} --filename ditch_seq_multiqc_report

# FASTQ statistics
if command -v seqkit &> /dev/null; then
    seqkit stats -j ${SLURM_CPUS_PER_TASK} ${FASTQ_SOURCE}/*.fastq.gz > ${SCRATCH}/processed_fastq_stats.txt
    seqkit stats -j ${SLURM_CPUS_PER_TASK} ${SCRATCH}/trimmed_fastq/*.fastq.gz > ${SCRATCH}/trimmed_fastq_stats.txt
fi

# BAM statistics
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

# Copy the log file to OUTPUT_BASE
cp /group/ll005/lditchburn/SLURM_Logs/slurm_log-${SLURM_JOB_NAME}-${SLURM_JOBID}.out ${OUTPUT_BASE}/

# Optionally clean up scratch
rm -rf $SCRATCH

################################################################################
# COMPLETION
################################################################################

echo ""
echo "============================================================"
echo "DITCH-seq bulk pipeline completed at $(date)"
echo "Results saved to: ${OUTPUT_BASE}"
echo "============================================================"
