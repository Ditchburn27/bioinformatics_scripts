for i in $(ls computeMatrix_output/*)
do
    echo $i
    output_file=$(basename ${i/.computeMatrix}_heatmap.png)
    plotHeatmap -m $i \
    -o plots/${output_file} \
    --colorList white,red \
    --sortRegions descend \
    --legendLocation upper-left \
    --whatToShow 'heatmap and colorbar'
done
