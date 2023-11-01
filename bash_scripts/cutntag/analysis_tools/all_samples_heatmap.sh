if [ -z "$1" ]
then
      echo -ne "
      ################\n
      First argument is empty. Please provide path to sample name text file.
      ################\n"
      exit 1
else
      sample_names=$1


      echo -ne "\nPath to sample name text file is\n\t${sample_names}\n"
fi

names=$(<$sample_names)

plotHeatmap -m computeMatrix_output/all_samples_computeMatrix.gz \
-o plots/all_samples_heatmap \
--colorList white,red \
--sortRegions descend \
--legendLocation upper-left \
--whatToShow 'heatmap and colorbar' \
--samplesLabel $names