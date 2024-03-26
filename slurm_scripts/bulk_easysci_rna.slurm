#!/bin/bash --login

#SBATCH --job-name=Bulk_easysci
#SBATCH --partition=peb
#SBATCH --mem-per-cpu=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --time=24:00:00   
#SBATCH --export=NONE
##SBATCH --mail-user=22720224@student.uwa.edu.au
##SBATCH --mail-type=BEGIN,END

# Start of job
echo $SLURM_JOB_NAME job started at  `date`

# To compile with the GNU toolchain
module load Anaconda3/2020.11
conda activate /group/ll005/envs/easy_sci

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
JOBNAME=${SLURM_JOB_NAME}

#  Take command line inputs
if [ -z "$1" ]
then
      echo -ne "
      ################\n
      First argument is empty. Please provide fastq_folder.\n
      ################\n"
      exit 1
else
      fastq_folder=$1


      echo -ne "\nPath to  path to fastq_folder is\n\t${fastq_folder}\n"
fi

if [ -z "$2" ]
then
      echo -ne "
      ################\n
      Second argument is empty. Please provide sample_ID file.\n
      ################\n"
      exit 1
else
      sample_ID=$2

fi

if [ -z "$3" ]
then
      echo -ne "
      ################\n
      Third argument is empty. Please provide output_folder name.\n
      ################\n"
      exit 1
else
      export output_folder=$3


      echo -ne "\nPath to  path to output_folder is\n\t${output_folder}\n"
fi

if [ -z "$4" ]
then
      echo -ne "
      ################\n
      Fourth argument is empty. Please provide star index.\n
      ################\n"
      exit 1
else
      index=$4


      echo -ne "\nPath to  path to star index is\n\t${index}\n"
fi

if [ -z "$5" ]
then
      echo -ne "
      ################\n
      Fith argument is empty. Please provide gtf_file.\n
      ################\n"
      exit 1
else
      gtf_file=$5


      echo -ne "\nPath to  path to gtf_file is\n\t${gtf_file}\n"
fi

if [ -z "$6" ]
then
      echo -ne "
      ################\n
      Sixth argument is empty. Please provide gtf_file_exons.\n
      ################\n"
      exit 1
else
      gtf_file_exons=$6


      echo -ne "\nPath to  path to gtf_file_exons is\n\t${gtf_file_exons}\n"
fi

if [ -z "$7" ]
then
      echo -ne "
      ################\n
      Seventh argument is empty. Please provide RT_barcode_matching_file.\n
      ################\n"
      exit 1
else
      RT_barcode_matching_file=$7


      echo -ne "\nPath to  path to RT_barcode_matching_file is\n\t${RT_barcode_matching_file}\n"
fi

if [ -z "$8" ]
then
      echo -ne "
      ################\n
      Eighth argument is empty. Please provide number of cores.\n
      ################\n"
      exit 1
else
      cores=$8


      echo -ne "\nNumber of cores is\n\t${cores}\n"
fi

## Barcode information files
script_folder=$(dirname "$0")/script_folder/
## define the location of the ligation barcodes
ligation_barcode_file=$script_folder/barcode_files/ligation_barcodes.pickle2
## define the location of the RT barcodes
RT_barcode_file=$script_folder/barcode_files/RT_barcodes.pickle2
## define the location of the randomN RT barcodes
randomN_barcode_file=$script_folder/barcode_files/RT_randomN_barcodes.txt

## Change the name of the fastq files to standard format
echo "Changing the name of the fastq files..."
for sample in $(cat $sample_ID); do echo changing name $sample; mv $fastq_folder/*$sample*R1*.fastq.gz $fastq_folder/$sample.R1.fastq.gz; \
mv $fastq_folder/*$sample*R2*.fastq.gz $fastq_folder/$sample.R2.fastq.gz; mv $fastq_folder/*$sample*R3*.fastq.gz $fastq_folder/$sample.R3.fastq.gz; done
echo

## Barcode the reads
echo "Barcoding reads..."
script=$script_folder/barcoding_reads_paired.py
mkdir -p $output_folder/barcoded_fastqs
python $script $fastq_folder $sample_ID $output_folder/barcoded_fastqs $ligation_barcode_file $RT_barcode_file $cores $randomN_barcode_file
echo "Done barcoding reads"

## Trim the reads
echo "Start trimming..."
mkdir -p $output_folder/trimmed_fastqs/

trim() {
    echo "Trimming sample: $1"
    trim_galore --paired $output_folder/barcoded_fastqs/$1*R1*.gz \
    $output_folder/barcoded_fastqs/$1*R2*.gz -a2 AAAAAAAA --stringency 3 -o $output_folder/trimmed_fastqs/
}

export -f trim

parallel -j $cores trim ::: $(cat $sample_ID)
echo "Done trimming..."
echo

## STAR alignment
echo "Start alignment using STAR..."
mkdir -p $output_folder/STAR_alignment
STAR --genomeDir $index --genomeLoad Remove

for sample in $(cat $sample_ID)
do 
echo "Aligning $sample..."
STAR --runThreadN $cores --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat \
 --readFilesIn $output_folder/trimmed_fastqs/$sample*R1*gz $output_folder/trimmed_fastqs/$sample*R2*gz \
 --outFileNamePrefix $output_folder/STAR_alignment/$sample --genomeLoad LoadAndKeep
done

STAR --genomeDir $index --genomeLoad Remove
echo "Done aligning"
echo


## Transfer STAR log files
mkdir -p $output_folder/report/Log_files
mv Log.out $output_folder/report/Log_files/
mv Log.progress.out $output_folder/report/Log_files/
mv Aligned.out.sam $output_folder/report/Log_files/
mv SJ.out.tab $output_folder/report/Log_files/
mv Log.final.out $output_folder/report/Log_files/


## Sorting and filtering SAM files
echo "Start sorting and filtering..."
mkdir -p $output_folder/filtered_sam

sort_filter() {
    echo "Sorting and filtering $1"
    samtools view -q 30 -f 2 -F 780 $output_folder/STAR_alignment/$1*.sam | sort -k1,1 -k3,3 -k4,4n > $output_folder/filtered_sam/$1.noheader.sam
    grep "@" $output_folder/STAR_alignment/$1*.sam > $output_folder/filtered_sam/$1.header.sam
    cat $output_folder/filtered_sam/$1.header.sam  $output_folder/filtered_sam/$1.noheader.sam > $output_folder/filtered_sam/$1.sam
    rm $output_folder/filtered_sam/$1.header.sam
    rm $output_folder/filtered_sam/$1.noheader.sam
    echo Filtering $1 done.
}
export -f sort_filter

parallel -j $cores sort_filter ::: $(cat $sample_ID)
echo "Done sorting and filtering"
echo


