#!/bin/bash

# script to move filtered feature matrix, to my directory and add library number prefix

for i in "$@"
do  cp /group/ll005/cellranger_out/$i/outs/per_sample_outs/$i/count/sample_filtered_feature_bc_matrix.h5 $MYGROUP/ | 
mv $MYGROUP/sample_filtered_feature_bc_matrix.h5 $MYGROUP/$i"_sample_filtered_feature_bc_matrix.h5"
done