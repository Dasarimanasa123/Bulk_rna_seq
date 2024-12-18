#!/bin/bash

# Define input and output directories
input_dir="deedup"
output_dir="BAM_files"

# Create output directory if it doesn't exist


# Loop over each FASTQ file in the input directory
for input_file in "$input_dir"/*.fastq.gz
do
    # Get the base filename without the directory path and extension
    base_name=$(basename "$input_file" .fastq.gz)

    # Define output file name
    output_file="$output_dir/${base_name}.bam"

    # Run fastp for each file
    
    hisat2 -x ref_data/mus_musculus_index -U $input_file  | samtools sort -o $output_file


    echo "Finished bamfile $input_file -> $output_file"

done

echo "All files processed."

