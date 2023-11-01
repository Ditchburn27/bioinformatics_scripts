computeMatrix reference-point \
--referencePoint center \
-S bigwig/RL*.sorted.rmDup.bw \
-R  MACS2/broad/RL*.sorted.rmDup.sorted_name_peaks.broadPeak \
-b 3000 -a 10000 \
-bs 50 \
-o computeMatrix_output/all_samples_computeMatrix.gz \
--missingDataAsZero \
--sortRegions descend \
-p 24
