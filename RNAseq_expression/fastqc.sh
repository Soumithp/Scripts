#!/bin/bash
module load fastqc/0.11.8
# Set the path to the folder containing fastq.gz files
input_folder="/work/SCCC/s226953/CGMH"

# Set the path to the folder where you want to store FastQC results
output_folder="/work/SCCC/s226953/CGMH/FastQC_results"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Run FastQC on all fastq.gz files in the input folder
for fastq_file in "$input_folder"/*.fastq.gz; do
    # Extract the file name without extension
    base_name=$(basename "$fastq_file" .fastq.gz)

    # Run FastQC
    fastqc "$fastq_file" -o "$output_folder"

    echo "FastQC completed for $base_name"
done
