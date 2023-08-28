#!/bin/bash

### Maps a plasmid against bowtie2 reference as well as de novo assembles fastq file
### requires: fastp, bowtie2, sambamba, SPAdes, samtools, seqkit, pigz, bcftools, tabix, unicycler, pilon, java
### version 0.7

version="0.7"
skipping="TRUE"
skipSpadesDenovoAssembly=0
qval=20
cleanSPADES=1
minLengthSPADES=0
uniqueID=""
runUnicycler=0

### Check dependencies

if ! [ -x "$(command -v fastp)" ]; then
  echo 'Error: fastp is not installed. Try installing with "conda install -c bioconda fastp' >&2
  exit 1
fi

if ! [ -x "$(command -v bowtie2-align-s)" ]; then
  echo 'Error: bowtie2-align-s is not installed. Try installing with "conda install -c bioconda bowtie2' >&2
  exit 1
fi

if ! [ -x "$(command -v sambamba)" ]; then
  echo 'Error: sambamba is not installed. Try installing with "conda install -c bioconda sambamba' >&2
  exit 1
fi
if ! [ -x "$(command -v pigz)" ]; then
  echo 'Error: pigz is not installed. Try installing with "conda install pigz' >&2
  exit 1
fi
if ! [ -x "$(command -v bcftools)" ]; then
  echo 'Error: bcftools is not installed. It should be part of samtools. Try installing with "conda install -c bioconda samtools' >&2
  exit 1
fi
if ! [ -x "$(command -v tabix)" ]; then
  echo 'Error: tabix is not installed. Try installing with "conda install -c bioconda tabix' >&2
  exit 1
fi
if ! [ -x "$(command -v plasmidspades.py)" ]; then
  echo 'Error: SPADES is not installed. Try installing with "conda install -c bioconda spades' >&2
  exit 1
fi
if ! [ -x "$(command -v seqkit)" ]; then
  echo 'Error: seqkit is not installed. Try installing with "conda install -c bioconda seqkit' >&2
  exit 1
fi
if ! [ -x "$(command -v unicycler)" ]; then
  echo 'Error: unicycler is not installed. Try installing with "conda install -c bioconda unicycler' >&2
  exit 1
fi
if ! [ -x "$(command -v pilon)" ]; then
  echo 'Error: pilon.jar is not installed. Try installing with "conda install -c bioconda pilon' >&2
  exit 1
fi




usage () { 
	cat << EOF
########## NEED SOME HELP ######## 	
plasmid_mapper_PE.sh ${version}
Please do not use special character for the names or spaces, "_" is acceptable.
OPTIONS:
    -h    Show this message
    -1    <file_r1.fastq> Fastq file, paired end read 1, gzipped is ok.
    -2    <file_r2.fastq> Fastq file, paired end read 2, gzipped is ok.
    -r    <ref.fa> Fasta reference file for mapping with bowtie2. Optional. Must end in *.fa !
    -s    skip SPADES plasmid de-novo assembly, default = FALSE
    -q    quality score for trimming, default = 20
    -c    clean up SPADES assembly directory and only save simplified_contigs that contain sequences, default = FALSE	
    -m    minimum length of SPADES assembly to be filtered out, default =  Disabled
    -u    unique ID to prepend file names with and fasta entries of SPADES nodes
    -y    unicycler assembly on trimmed fastq files

Example: 
bash plasmid_mapper_PE.sh \\
		-1 RL1093_2018_04_18_INV_pLenti_dCas9_SunTag_2_5059_S6_R1.fastq \\
		-2 RL1093_2018_04_18_INV_pLenti_dCas9_SunTag_2_5059_S6_R2.fastq \\
		-r INV_pLenti-tetO_NheI_HA-dCas9-SunTag_MluI_IRES_BFP_bGHpolyA_PacI_pEf1a_Puro-T2A-tTA-VP16_GG_ligation.fa \\
		-c
		-m 300
		-u RL1093
		-y

Script written by Christian Pflueger
EOF
	exit
}

#parse the options
while getopts ":1:2:r:sycq:m:u:" o; do
	case "${o}" in
		1)
			fastQ_f=${OPTARG}
			;;
		2)
			fastQ_r=${OPTARG}
			;;
		r)
        	fastaRef=${OPTARG}
			;;
		s)
			skipSpadesDenovoAssembly=1
			;;
		c)
			cleanSPADES=0
			;;
		q)
			qval=${OPTARG}
			;;
		m)
			minLengthSPADES=${OPTARG}
			;;
		u)
			uniqueID=${OPTARG}
			;;
		y)
			runUnicycler=1
			;;
		*)
			usage
            ;;
    esac
done
shift $((OPTIND-1))


#check the options
if [ -z "${fastQ_f}" ]; then
    echo
    echo "Read 1 fastq file is missing"
    usage
fi
if [ -z "${fastQ_r}" ]; then
    echo
    echo "Read 2 fastq file is missing"
    usage
fi
if [ -z "${fastaRef}" ]; then
    echo
    echo "Path/to/fasta.fa file is missing. Skipping reference assembly..."
    skipping="TRUE"
elif [ "${fastaRef}" == "na" ]; then
	# no reference availabe -> specifically set by user with na -> skip reference assembly
	echo
	echo "No fasta reference file given by providing option '-r na'. Skipping reference assembly..."
	echo
	skipping="TRUE"
else
	if [ "${fastaRef##*.}" != "fa" ] && [ "${fastaRef##*.}" != "fasta" ]; then
    	echo
    	echo "####################################################"
    	echo "Fasta file does not end in either '.fa' or '.fasta' "
    	echo "####################################################"
    	echo
    	usage
	fi
	skipping="FALSE"
fi

fastaRefFile=${fastaRef}
bowtieIndex=${fastaRefFile%.fa*}

# 1. Trim Reads (Nextera Adapters)

fastp -w 4 -q ${qval} \
	  -i ${fastQ_f} \
	  -I ${fastQ_r} \
	  -o "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" \
	  -O "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq" \
	  -h "${fastQ_f%_S[0-9]*}_fastp_report.html" \
	  -j "${fastQ_f%_S[0-9]*}_fastp_report.json" \
	  -R "${fastQ_f%_S[0-9]*}_fastp_report" \
	  --adapter_sequence CTGTCTCTTATACAC \
	  --adapter_sequence_r2 TGTATAAGAGACAG \
	  -x -y

### Reference mapping if fasta reference file was provided

if [ "${skipping}" == "FALSE" ]; then
	# 2. Build Bowtie2 Index

	bowtie2-build -f ${fastaRef} ${bowtieIndex} &> "${fastQ_f%_S[0-9]*}_bt2_index.log"

	# 3. Bowtie2 map against index

	bowtie2-align-s -x ${bowtieIndex} \
					--very-sensitive  \
					-1 "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" \
					-2 "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq" \
					-S "${fastQ_f%_S[0-9]*}.sam" &> "${fastQ_f%_S[0-9]*}_bowtie2.log"
	echo
	cat "${fastQ_f%_S[0-9]*}_bowtie2.log"
				
	# 4. Sam to Bam conversion, filtering and removal of multimapper

	sambamba view -q -F "not(secondary_alignment)" -S -f bam "${fastQ_f%_S[0-9]*}.sam" |\
	sambamba sort -q /dev/stdin -o "${fastQ_f%_S[0-9]*}.bam"

fi

# 5. Run de-novo assembly with plasmids_spades

if [ ${skipSpadesDenovoAssembly} == 0 ]; then
	
	# create output directory
	spades_dir="${fastQ_f%_S[0-9]*}_SPADES_assembly"
	mkdir ${spades_dir}

	# run plasmid spades for de-novo assembly
	plasmidspades.py \
	-1 "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" \
	-2 "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq" \
	-o ${spades_dir} &> "${spades_dir}/${spades_dir}.log"

	# remove temp directory
	rm -r "${spades_dir}/tmp"

	# SPADES cleanup option

	if [ ${cleanSPADES} == 0 ]; then

		mkdir ${spades_dir/_assembly/_essential}

		# add the Kmer info to simplified fasta file name
		find . -name "simplified_contigs.fasta" -size +1 | while read line; do \
			n=${line/\/simplified/_simplified}; \
			n=${n##*/}; \
			mv ${line} ${line%/*}/${n}; \
		done

		# copy log files and simplified fasta SPADES assemblies to 'essential' folder
		find ${spades_dir} \( -name "*simplified_contigs.fasta" -o -name "*.log" \) \
		-exec cp {} ${spades_dir/_assembly/_essential} \;

		# remove the SPADES_assembly folder
		rm -r ${spades_dir}

		# Filter nodes for minimum length assemblies 
		if [ ${minLengthSPADES} > 0 ]; then

			# gather simplified SPADES fasta files
			for f in ${spades_dir/_assembly/_essential}/*simplified_contigs.fasta; do
				
				# output file name

				outName="${f%_simplified_contigs.*}_SPADES_${minLengthSPADES}bp.fasta"

				# filter length of fasta files
				seqkit seq -m ${minLengthSPADES} ${f} > ${outName}
			
				# check if unique ID was set by user
				if [[ ! -z "$uniqueID" ]]; then

					# prepend unique identifier variable to file name	
					outNameUnique=${outName/\/K/\/${uniqueID}_K}

					# rename file name with prepended uniqueID
					mv ${outName} ${outNameUnique} #all plasmid spades files start with Kxx for the kmer number

					# rename fasta node entries
					export uniqueID
					perl -i -pe 's/\>/\>$ENV{uniqueID}_/g' ${outNameUnique}
				fi

			done

		fi
	fi

fi

# 6. Unicycler de novo assembly
if [ "${runUnicycler}" == 1 ]; then

	# set output directory
	unicycler_dir="${fastQ_f%_S[0-9]*}_unicycler_assembly"
	mkdir ${unicycler_dir}

	# run unicyle for de-novo assembly
	unicycler \
	-1 "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" \
	-2 "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq" \
	-t 4 \
	--mode normal \
	--keep 0 \
	--verbosity 0 \
	-o ${unicycler_dir} &> "${unicycler_dir}/${unicycler_dir}.log"

	# rename "assembly.fasta" to a meaningful file name
	mv "${unicycler_dir}/assembly.fasta" "${unicycler_dir}/${fastQ_f%_S[0-9]*}_unicycler_assembly.fasta"
	mv "${unicycler_dir}/assembly.gfa" "${unicycler_dir}/${fastQ_f%_S[0-9]*}_unicycler_assembly.gfa"
	mv "${unicycler_dir}/unicycler.log" "${unicycler_dir}/${unicycler_dir}.log"

fi



# 7. Create consensus fasta
if [ "${skipping}" == "FALSE" ]; then
	
	bcftools mpileup -Ou -f ${fastaRef} "${fastQ_f%_S[0-9]*}.bam" |\
	bcftools call -mv -Oz -o "${fastaRef/.fa/.vcf.gz}"
	tabix "${fastaRef/.fa/.vcf.gz}"
	cat ${fastaRef} | bcftools consensus "${fastaRef/.fa/.vcf.gz}" > "${fastaRef/.fa/_consensus.fa}"

	# check if unique ID was set by user
	if [[ ! -z "$uniqueID" ]]; then

		# prepend unique identifier variable to consensus fasta
		mv ${fastaRef/.fa/_consensus.fa} "${uniqueID}_${fastaRef/.fa/_consensus.fa}"

	fi
fi


# 8. Assemble unmapped reads with unicycler to discover possible plasmid mixtures or other contaminants
# Only works if fasta reference was provided

if [ "${skipping}" == "FALSE" ]; then

	# create directory for unmapped files
	unmapped_reads_dir="${fastQ_f%_S[0-9]*}_unmapped_reads"
	mkdir ${unmapped_reads_dir}

	# extract unmapped reads from bam file
	samtools view -@ 2 -f4 -bh "${fastQ_f%_S[0-9]*}.bam" > "${fastQ_f%_S[0-9]*}.unmapped.bam"

	# convert unmapped bam to fastq files
	samtools fastq -@ 2 "${fastQ_f%_S[0-9]*}.unmapped.bam" \
	-1 "${fastQ_f%_S[0-9]*}.unmapped.R1.fastq.gz" \
	-2 "${fastQ_f%_S[0-9]*}.unmapped.R2.fastq.gz" \
	-0 /dev/null -s /dev/null -n

	# assemble the unmapped reads with unicycler
	unicycler \
	-1 "${fastQ_f%_S[0-9]*}.unmapped.R1.fastq.gz" \
	-2 "${fastQ_f%_S[0-9]*}.unmapped.R2.fastq.gz" \
	-t 4 \
	--mode normal \
	--keep 0 \
	--verbosity 0 \
	-o ${unmapped_reads_dir} &> "${unmapped_reads_dir}/${unmapped_reads_dir}.log"

	# rename "assembly.fasta" to a meaningful file name
	mv "${unmapped_reads_dir}/assembly.fasta" "${unmapped_reads_dir}/${fastQ_f%_S[0-9]*}_unicycler_assembly.fasta"
	mv "${unmapped_reads_dir}/assembly.gfa" "${unmapped_reads_dir}/${fastQ_f%_S[0-9]*}_unicycler_assembly.gfa"
	mv "${unmapped_reads_dir}/unicycler.log" "${unmapped_reads_dir}/${unicycler_dir}.log"

	# cleanup
	mv "${fastQ_f%_S[0-9]*}.unmapped.R1.fastq.gz" \
	"${fastQ_f%_S[0-9]*}.unmapped.R2.fastq.gz" \
	"${fastQ_f%_S[0-9]*}.unmapped.bam" \
	${unmapped_reads_dir}


fi

# 9. Cleanup

pigz -p 12 "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq"
find . -name "*.sam" -delete
find . -name "*.bt2" -delete
find . -name "0" -delete