#!/bin/bash

# Define input and output directories
input_dir="trimm_data"
output_dir="deedup"

# Create output directory if it doesn't exist


# Loop over each FASTQ file in the input directory
for input_file in "$input_dir"/*.fastq.gz
do
    # Get the base filename without the directory path and extension
    base_name=$(basename "$input_file" .fastq.gz)

    # Define output file name
    output_file="$output_dir/${base_name}_deedup.fastq.gz"

    # Run fastp for each file
    
    fastp -i "$input_file" -o "$output_file" --dedup --thread 4


    echo "Finished deedup $input_file -> $output_file"

done

echo "All files processed."

