#!/bin/bash --login

#SBATCH --job-name=CUTnTag
#SBATCH --partition=peb
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --time=8:00:00   
#SBATCH --export=NONE
##SBATCH --mail-user=22720224@student.uwa.edu.au
##SBATCH --mail-type=BEGIN,END

set -eu -o pipefail -o verbose

# Start of job
echo $SLURM_JOB_NAME job started at  `date`

# To compile with the GNU toolchain
module load Anaconda3/2020.11
conda activate /group/ll005/envs/cutntag

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
JOBNAME=${SLURM_JOB_NAME}
SCRATCH=$MYSCRATCH/$JOBNAME/$SLURM_JOBID
BOWTIE2_INDEX=/group/ll005/reference/bowtie2_hg38/hg38
BLACKLIST=/group/ll005/reference/bowtie2_hg38/blacklist/hg38-blacklist.v2.bed
###############################################
# Creates a unique directory in the SCRATCH directory for this job to run in.
if [ ! -d $SCRATCH ]; then 
    mkdir -p $SCRATCH 
fi 
echo SCRATCH is $SCRATCH

#############################################
#   Copy input files to $SCRATCH
#   then change directory to $SCRATCH

mkdir -p ${SCRATCH}/{fastq,fastqc,trimmed_fastq,sorted_bam,dedup_bam,clean_bed,bigwig,plots,fragment_bed}
mkdir -p ${SCRATCH}/MACS2/{narrow,broad}
mkdir -p ${SCRATCH}/logs/{fastp,bowtie2,dedup}

### COPY files to scratch
cp fastq/*gz ${SCRATCH}/fastq

## IMPORTANT: change directory to scratch

cd $SCRATCH
###################
#1. FastqQC

fastqc fastq/*fastq.gz -o ${SCRATCH}/fastqc --threads ${SLURM_CPUS_PER_TASK}

#2.get the read1 and read2 name and put in a file
ls fastq/*_R1*fastq.gz >> read1
ls fastq/*_R2*fastq.gz >> read2

#3.Adapter Trimming with fastp

for i in $(cat read1)
do 
	read j
	echo $i $j
	fastp_R1=$(basename ${i/.fastq.gz}_trimmed.fastq.gz)
	fastp_R2=$(basename ${j/.fastq.gz}_trimmed.fastq.gz)
	html=$(basename ${i%_R1_001.fastq.gz}_fastp_report.html)
	report=$(basename ${i%_R1_001.fastq.gz}_fastp_report)
	fastp_log=$(basename ${i%_R1_001.fastq.gz}_fastp.out)
	fastp -w ${SLURM_CPUS_PER_TASK} -x -y \
	-i $i \
	-I $j \
	-o ${SCRATCH}/trimmed_fastq/${fastp_R1} \
	-O ${SCRATCH}/trimmed_fastq/${fastp_R2} \
	-h ${SCRATCH}/trimmed_fastq/${html} \
	-R ${SCRATCH}/trimmed_fastq/${report} &> ${SCRATCH}/logs/fastp/${fastp_log}
done < read2

#4.get the read1 and read2 name of trimmed fastq and put in a file

ls trimmed_fastq/*_R1*fastq.gz >> trimmed_read1
ls trimmed_fastq/*_R2*fastq.gz >> trimmed_read2

#5.mapping with bowtie2 removing unmapped reads and convert to sorted bam and create index
for i in $(cat trimmed_read1)
do 
	read j
	echo $i $j
	bowtie2_log=$(basename ${i%_R1_001_trimmed.fastq.gz}_bowtie2.txt)
	sorted_bam=$(basename ${i%_R1_001_trimmed.fastq.gz}.sorted.bam)
	bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 25 -x ${BOWTIE2_INDEX} -1 $i -2 $j 2> logs/bowtie2/${bowtie2_log} | samtools view -bS -F 4 - | samtools sort -o sorted_bam/${sorted_bam} - 
done < trimmed_read2
parallel -j5 "samtools index {} {.}.bam.bai" ::: sorted_bam/*bam

#6.mark duplicates with picard tools

parallel -j15 picard MarkDuplicates -I {} -O dedup_bam/{/.}.rmDup.bam --REMOVE_DUPLICATES true -M logs/dedup/{/.}_picard.rmDup.txt ::: $(ls sorted_bam/*bam)
parallel -j5 "samtools index {} {.}.bam.bai" ::: dedup_bam/*.sorted.rmDup.bam

#7. Sort BAM files by name using parallel processing
parallel --jobs 10 "samtools sort -@ 3 -n {} -o fragment_bed/{/.}.sorted_name.bam" ::: $(ls dedup_bam/*.bam)


#8.Convert the sorted BAM file to fragment files 
# adjusting Tn5 cut by moved forward by 4 bp from a left-most alignment position and backward 5 bp from the right-most alignment position'
# remove chrM
# remove fragments overlapping with blacklist regions

parallel -j10 "bedtools bamtobed -bedpe -i {} | cut -f 1,2,6 | awk -v OFS="\t" '{print $1;$2+4;$3-5}'| grep -v "chrM" | bedtools intersect -v -a - -b $BLACKLIST > fragment_bed/{/.}.bed" ::: $(ls fragment_bed/*.bam)
rm fragment_bed/*.bam

#9.Call peaks from fragments files

parallel -j5 macs2 callpeak -t {} -f BEDPE --keep-dup all --broad-cutoff 0.05 --broad --gsize hs -n {/.} --outdir MACS2/broad ::: $(ls fragment_bed/*.bed)
parallel -j5 macs2 callpeak -t {} -f BEDPE --keep-dup all -q 0.05 --gsize hs -n {/.} --outdir MACS2/narrow ::: $(ls fragment_bed/*.bed)


#10.Create bigwig files

parallel -j2 bamCoverage -p 10 -b {} --normalizeUsing CPM --blackListFileName $BLACKLIST -o bigwig/{/.}.bw ::: $(ls dedup_bam/*.bam)


#11.Multiqc
multiqc .

#12.summary of fastq
seqkit -j 10 stats fastq/*.fastq.gz > fastq/raw_fastq_summary.txt
seqkit -j 10 stats trimmed_fastq/*.fastq.gz > trimmed_fastq/trimmed_fastq_summary.txt

#12.Summary of bam files
for i in sorted_bam/*bam;do echo -ne "$i\t" >> bam_stats_summary.txt; samtools view -@ 10 -F 0x40 $i | cut -f1 | sort | uniq | wc -l >> bam_stats_summary.txt;echo;done 



#############################################
# Move all files from scratch to data folders
#rsync -avzc * 

###########################
# Clean up $SCRATCH 

# Save files that were created
#mv ${SCRATCH}/* ${RESULTS}

# Delete SRATCH
#rm -r $SCRATCH

echo $JOBNAME job finished at  `date`