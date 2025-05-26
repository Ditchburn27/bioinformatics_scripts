#!/bin/bash

# Script for filtering out redundant exons from gtf file. 
# Used for easy-sci RNA-seq pipeline for exon counting step. 

gtf_file=$1
output=$2

awk '
{
    if ($0 ~ /exon_id/) {
        match($0, /exon_id "([^"]+)"/, arr)
        exon_id = arr[1]
        if (!(exon_id in seen)) {
            seen[exon_id] = 1
            print $0
        }
    } else {
        print $0
    }
}' $gtf_file > $output