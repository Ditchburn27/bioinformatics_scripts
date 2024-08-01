#!/bin/bash

# Input and output file names
input_file=$1
output_file=$2

# Clear the output file if it already exists
> "$output_file"

# Process each line of the input file
while IFS=$'\t' read -r -a fields; do
    # Get the last field (attribute field)
    attributes=${fields[8]}
    
    # Remove the trailing semicolon and add spaces after semicolons
    formatted_attributes=$(echo "$attributes" | sed 's/; /;/g' | sed 's/;/; /g' | sed 's/; $//')

    # Reconstruct the line with formatted attributes
    fields[8]=$formatted_attributes
    formatted_line=$(printf "%s\t" "${fields[@]}")
    formatted_line=${formatted_line%$'\t'}
    
    # Write the formatted line to the output file
    echo -e "$formatted_line" >> "$output_file"
done < "$input_file"