
setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/07.23.2024_Epacad_rat_analysis")
# Load required library
library(dplyr)

# Read input data
df <- read.table("CDAHFD_Epacad_rat_RLE.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Count genes with non-zero values in each sample (column)
gene_counts <- colSums(df > 0)

# Print or write gene counts (modify as needed)
print(gene_counts)
write.table(gene_counts, file = "output_epacad_rat_gene_count.txt", sep = "\t", quote = FALSE)
