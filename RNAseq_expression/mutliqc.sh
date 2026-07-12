#!/bin/bash
# Set the path to your main folder containing FastQC files
main_folder="/work/SCCC/soumith/Epacad_rat"

# Set the output directory for MultiQC results
output_dir="/work/SCCC/soumith/Epacad_rat"
# Set the name of the Excel file to store MultiQC results
excel_file="multiqc_results.xlsx"
multiqc --dirs -o "$output_dir" "$main_folder"
