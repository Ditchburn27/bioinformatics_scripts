#!/bin/bash

# Check if dir_list.txt exists
if [ ! -f "dir_list.txt" ]; then
    echo "Error: dir_list.txt not found."
    exit 1
fi

# Function to process fragments for a directory
process_directory() {
    directory="$1"
    cd "$directory" || return
    
    # Navigate to the outs folder
    outs_folder="outs"
    if [ ! -d "$outs_folder" ]; then
        echo "Error: $outs_folder folder not found in $directory."
        cd ..
        return
    fi
    cd "$outs_folder" || return
    
    # Execute the prefix.sh script
    if [ -f "fragments.tsv.gz" ]; then
        echo "Processing fragments in $directory/outs..."
        /group/ll005/processed_data/cell_braindev_paper_mapped_snATAC/prefix.sh "fragments.tsv.gz"
    else
        echo "Error: fragments.tsv.gz not found in $directory/outs."
    fi
    
    # Navigate back to the original directory
    cd ../..
}

# Export function for GNU Parallel
export -f process_directory

# Process directories in parallel
cat dir_list.txt | parallel process_directory
