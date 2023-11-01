mkdir computeMatrix_output
ls bigwig/* >> bigwig_files
ls MACS2/broad/RL*.sorted.rmDup.sorted_name_peaks.broadPeak >> peak_files

for i in $(cat bigwig_files)
do
    read j
    echo $i $j 
    output_file=$(basename ${i/.sorted.rmDup.bw}.computeMatrix)
    computeMatrix reference-point \
    --referencePoint center \
    -S $i \
    -R $j \
    -b 3000 -a 10000 \
    -bs 50 \
    -o computeMatrix_output/${output_file} \
    --missingDataAsZero \
    --sortRegions descend \
    -p 24
done < peak_files
