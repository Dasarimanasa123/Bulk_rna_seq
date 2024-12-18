#!/bin/bash

# Define input and output directories
input_dir="raw_data"
output_dir="trimm_data"

# Create output directory if it doesn't exist


# Loop over each FASTQ file in the input directory
for input_file in "$input_dir"/*.fastq.gz
do
    # Get the base filename without the directory path and extension
    base_name=$(basename "$input_file" .fastq.gz)

    # Define output file name
    output_file="$output_dir/${base_name}_trimmed.fastq.gz"

    # Run Trimmomatic for each file
    trimmomatic SE -threads 4 "$input_file" "$output_file" \
    HEADCROP:5 LEADING:20 MINLEN:50 SLIDINGWINDOW:4:15

    echo "Finished trimming $input_file -> $output_file"

done

echo "All files processed."

