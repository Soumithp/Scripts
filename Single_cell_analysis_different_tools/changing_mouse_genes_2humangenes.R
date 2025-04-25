###to change mouse genes to human genes
setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/")
# Load necessary libraries
library(dplyr)
library(tidyr)

# Load the gene signatures table
gene_signatures <- read.delim("zonation_HSC.txt", header = TRUE, stringsAsFactors = FALSE)

# Load the reference file (Probe Set ID and Gene Symbol)
reference_table <- read.delim("Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v2023.2.Hs.txt", header = TRUE, stringsAsFactors = FALSE)

# Initialize an empty list to store the results
result <- list()

# Iterate over each column in the gene signature table
for (col in colnames(gene_signatures)) {
  
  # Extract the gene list from the column
  genes <- gene_signatures[[col]]
  
  # Initialize a vector to store the translated genes
  translated_genes <- c()
  
  # Iterate over each gene in the list
  for (gene in genes) {
    
    # Split gene names that are separated by ";"
    split_genes <- unlist(strsplit(gene, ";"))
    
    # Find corresponding gene symbols in the reference table using "Probe Set ID"
    translated_gene <- sapply(split_genes, function(g) {
      match_idx <- match(g, reference_table$Probe.Set.ID)  # Case-sensitive search in Probe Set ID
      if (!is.na(match_idx)) {
        return(reference_table$Gene.Symbol[match_idx])  # Return the Gene Symbol if a match is found
      } else {
        return(g)  # Return the original gene if no match is found
      }
    })
    
    # Combine the translated gene names (if any)
    translated_genes <- c(translated_genes, paste(translated_gene, collapse = ";"))
  }
  
  # Add the translated column to the result
  result[[col]] <- translated_genes
}

# Convert the result list into a data frame
final_result <- as.data.frame(result)

# Write the result to a new text file
write.table(final_result, "translated_gene_signatures_zonation_08.22.2024.txt", sep = "\t", row.names = FALSE, quote = FALSE)
