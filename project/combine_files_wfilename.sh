#!/bin/bash

# Specify the input directory and output file
input_directory=~/finalProject/myresults/per_contig_row5
output_file="combined_file.txt"

# Remove existing output file if it exists
[ -e "$output_file" ] && rm "$output_file"

# Iterate through all files in the input directory
for file in "/users/n/h/nheaphy/finalProject/myresults/per_contig_row5"/*; do
    # Get the filename without the path
    filename=$(basename "$file")

    # Append the filename and content as a new row to the output file
    awk -v filename="$filename" '{print filename "\t" $0}' "$file" >> "$output_file"
done

echo "Combined file created: $output_file"
