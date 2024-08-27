#! bin/bash
### Script to take a gtf file and extract only the exons lines.
### exon gtf file is required for easy-sci RNA-seq pipeline.
### The feature column is changed from 'exon' to 'gene'.
### The final column is reformatted to match the example
### provided by easy-sci (see their github) 
input_gtf=$1
output_gtf_name=$2

echo "Extracting exon entries from gtf..."
awk -F "\t" '$3=="exon"' $input_gtf > tmp1.gtf
tr ';' '\t' < tmp1.gtf > tmp2.gtf
awk -F'"' 'BEGIN {OFS="\t"} {$1=$1; gsub(/ +/, "\t"); gsub(/"/, "", $2); print}' tmp2.gtf > tmp3.gtf
awk -F'\t' 'BEGIN {OFS="\t"} {$9 = $9 "\t" substr($11, 1, length($11)) "-" substr($54, 1, length($54)); print}' tmp3.gtf > tmp4.gtf
cut -f 1-10,30-37 tmp4.gtf > tmp5.gtf
awk -v field=10 'BEGIN {OFS="\t"} {if (NF >= field) $field = "\"" $field "\""; print}' tmp5.gtf > tmp6.gtf
awk -v field=12 'BEGIN {OFS="\t"} {if (NF >= field) $field = "\"" $field "\""; print}' tmp6.gtf > tmp7.gtf
awk -v field=14 'BEGIN {OFS="\t"} {if (NF >= field) $field = "\"" $field "\""; print}' tmp7.gtf > tmp8.gtf
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9 " " $10 ";" $11 " " $12 ";"  $13 " " $14 ";" }' tmp8.gtf > tmp9.gtf
awk 'BEGIN {OFS="\t"} {if ($3 == "exon") $3 = "gene"; print}' tmp9.gtf > tmp10.gtf
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9 " " $10 " " $11 " " $12}' tmp10.gtf > $output_gtf_name

echo "exon gtf saved to $output_gtf_name."

# Clean up temporary files
rm tmp*.gtf