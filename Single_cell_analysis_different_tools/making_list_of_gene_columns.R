#making list from column 
setwd("C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki")
gene_signatures <- read.table("C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiverMP_sig_allcombined_wnodup.txt", header = TRUE, sep = "\t", check.names = FALSE)

# Specify the output file path
output_file <- "C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiver_MP_allsign_wodup_list_format.txt"

# Open the file connection
file_conn <- file(output_file, "w")

# Create a list to store each cell type with its genes as a comma-separated string
lapply(names(gene_signatures), function(cell_type) {
  # Extract genes for the cell type
  genes <- gene_signatures[[cell_type]]
  
  # Remove any NA values (if present) to clean up the list
  genes <- genes[!is.na(genes)]
  
  # Create a comma-separated string without quotation marks
  gene_string <- paste(genes, collapse = ", ")
  
  # Write the result with the cell type name to the file
  writeLines(paste0(cell_type, ": ", gene_string), con = file_conn)
  writeLines("\n", con = file_conn)  # Add a new line for better readability
})

# Close the file connection
close(file_conn)

cat("Gene signatures have been saved to:", output_file, "\n")