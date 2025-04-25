
library(Seurat)
library(tidyverse)
library(hdf5r)
library(ggplot2)
library(dplyr)
library(pheatmap)
setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/")


scLiver_MP_signature <- read.delim("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiverAtlas_allsign.txt")
counts <- Read10X_h5("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/hg38/hg38_filtered_feature_bc_matrix.h5")
df <- CreateSeuratObject(counts = counts)
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
df <- subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
df <- NormalizeData(df)
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(df)
df <- ScaleData(df, features = all.genes)
df <- RunPCA(df, features = VariableFeatures(object = df))
ElbowPlot(df)
#clustering
df <- FindNeighbors(df)
df <- FindClusters(df, resolution = 1)
df <- RunUMAP(df, dims = 1:10)
DimPlot(df, reduction = "umap", pt.size = 3)


signatures = list()

for(cell_type in colnames(scLiver_MP_signature)){
  print(paste('Extracting the cell type:',cell_type,'...'))
  signature_list = scLiver_MP_signature[[cell_type]]
  signatures[[cell_type]] = signature_list[which(signature_list!="")]
}
for(i in 1:ncol(scLiver_MP_signature)){
  print(i)
  df <- AddModuleScore(object = df, features = list(signatures[[i]]), name = paste0('Celltype_',names(signatures)[i]))
}

scores= df@meta.data[, grep("^Celltype_", colnames(df@meta.data))]
head(scores)
df$cell_type = apply(scores, 1, function(x) {
  cell_types <- names(signatures)
  cell_types[which.max(x)]
})


###cheking how many genes are present in the matrix and comparing with the gene signature
# Step 1: Extract the gene names from the expression matrix
genes_in_matrix <- rownames(df)
# Step 2: Initialize a vector to store the counts of genes present for each signature
genes_present_count <- c()
# Step 3: Loop through each signature to check how many genes are present in the matrix
for(cell_type in colnames(scLiver_MP_signature)){
  # Extract the list of genes for the current cell type
  signature_genes <- scLiver_MP_signature[[cell_type]]
  
  # Remove empty entries
  signature_genes <- signature_genes[which(signature_genes != "")]
  # Check which genes from the signature are present in the matrix
  present_genes <- signature_genes %in% genes_in_matrix
  # Count the number of present genes
  count_present <- sum(present_genes)
  # Store the count in the vector
  genes_present_count <- c(genes_present_count, count_present)
  # Print the result for each cell type
  print(paste('Cell type:', cell_type, '- Genes present:', count_present, '/', length(signature_genes)))
}

# Step 4: Create a data frame for visualization
gene_presence_df <- data.frame(
  Cell_Type = colnames(scLiver_MP_signature),
  Genes_Present = genes_present_count,
  Total_Genes = sapply(signatures, length)
)


#####using intersect
# Step 1: Calculate the intersection and store counts
gene_presence_df <- data.frame(
  Cell_Type = colnames(scLiver_MP_signature),
  Genes_Present = sapply(colnames(scLiver_MP_signature), function(cell_type) {
    length(intersect(scLiver_MP_signature[[cell_type]][scLiver_MP_signature[[cell_type]] != ""], genes_in_matrix))
  }),
  Total_Genes = sapply(signatures, length)
)

# Step 2: Plot the number of genes present in the matrix for each signature
ggplot(gene_presence_df, aes(x = Cell_Type, y = Genes_Present)) +
  geom_bar(stat = "identity", fill = "red", width = 0.4) +  # Reduce the bar width
  geom_text(aes(label = paste(Genes_Present, "/", Total_Genes)), vjust = -0.5, size = 3) +  # Reduce font size
  labs(title = "Presence of Signature Genes in Expression Matrix",
       x = "Cell Type",
       y = "Number of Genes Present") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),  # Rotate x-axis labels and reduce size
    axis.text.y = element_text(size = 5),  # Reduce y-axis label size
    plot.title = element_text(size = 10)   # Adjust title size
  )

