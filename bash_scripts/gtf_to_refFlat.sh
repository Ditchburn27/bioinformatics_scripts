#! /bin/bash
### Script to convert gtf files to refFlat format
### required by Picard CollectRnaSeqMetrics

gtf_file=$1
refFlat_file=$2

gtfToGenePred -genePredExt -geneNameAsName2 $gtf_file refFlat.tmp.txt
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > $refFlat_file
rm refFlat.tmp.txt