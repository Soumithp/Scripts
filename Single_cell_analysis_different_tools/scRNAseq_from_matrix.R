install.packages("Seurat")
install.packages("dplyr")
install.packages("fastmap")


library(Seurat)
library(tidyverse)
library(hdf5r)
library(ggplot2)
library(dplyr)
library(Matrix)
library(openxlsx)

setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/")
#reading the filtered data
svr_data<- Read10X_h5(filename= "C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/hg38/hg38_filtered_feature_bc_matrix.h5")
svr.seurat.obj<- CreateSeuratObject(counts= svr_data,  project = "Svr", min.cells= 3, min.features= 200)
##Qc %mt reads
svr.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(svr.seurat.obj, pattern = "^MT-")
VlnPlot(svr.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(svr.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
# 2. Filtering
svr.seurat.obj <- subset(svr.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
# 3. Normalize data 
svr.seurat.obj <- NormalizeData(svr.seurat.obj)

# # Step 4: Extract the normalized data matrix
# normalized_matrix <- GetAssayData(svr.seurat.obj, layer = "data")
# # Step 5: Save the matrix to a text file
# write.table(as.matrix(normalized_matrix), file = "ev59_filtered_normalized_matrix.txt", sep = "\t", quote = FALSE, col.names = NA)

# 4.Identify highly variable features 
svr.seurat.obj <- FindVariableFeatures(svr.seurat.obj, selection.method = "vst", nfeatures = 1000)

# 5. Scaling
all.genes <- rownames(svr.seurat.obj)
svr.seurat.obj <- ScaleData(svr.seurat.obj, features = all.genes)
# 6. Perform Linear dimensionality reduction
svr.seurat.obj <- RunPCA(svr.seurat.obj, features = VariableFeatures(object = svr.seurat.obj, npcs= 50))
# 7. Clustering
svr.seurat.obj <- FindNeighbors(svr.seurat.obj, dims = 1:20)
svr.seurat.obj <- FindClusters(svr.seurat.obj, resolution = 0.5)

# setting identity of clusters
Idents(svr.seurat.obj) <- "RNA_snn_res.0.5"
#Idents(veh.seurat.obj)
svr.seurat.obj <- RunUMAP(svr.seurat.obj, dims = 1:20)
DimPlot(svr.seurat.obj, reduction = "umap")


#############annotation

#########checking the expression 
########first to make the gene signatures into a list
# Step 1: Read and Prepare Gene Signatures
gene_signatures <- read.table("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/DDR1hepg2_FDR025_sig_20231018_FINAL.txt", header = TRUE, sep = "\t", check.names = FALSE)

# Convert gene signatures to a list of gene sets
gene_list <- lapply(gene_signatures, as.character)
features_list <- lapply(colnames(gene_signatures), function(x) gene_signatures[[x]])
names(features_list) <- colnames(gene_signatures)

rownames(svr.seurat.obj) <- toupper(rownames(svr.seurat.obj))

#########################identifying the expressed genes
# Step 2: Get the Count Matrix
counts_matrix <- GetAssayData(svr.seurat.obj, assay = "RNA", layer = "counts")

# Filter features_list to only include genes that are present in the counts_matrix
filtered_features_list <- lapply(features_list, function(genes) {
  intersect(genes, rownames(counts_matrix))
})

##########to get matreix of genes presetn or absent
# combined_list <- lapply(names(cells_with_signature), function(signature) {
#   df <- as.data.frame(cells_with_signature[[signature]])
#   df$GeneSignature <- signature
#   df
# })
# presence_absence_df <- do.call(rbind, combined_list)
# colnames(presence_absence_df) <- c(colnames(counts_matrix), "GeneSignature")
# 
# # Write to a text file
# write.table(presence_absence_df, "gene_signature_presence_absence.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# # Write to an Excel file
# write.xlsx(presence_absence_df, "gene_signature_presence_absence.xlsx")

####to get only gene names which are present in the cells
# Initialize list to store presence/absence data
# Combine all data frames into one

# Combine presence/absence matrices for each gene signature
combined_list <- lapply(names(cells_with_signature), function(signature) {
  genes <- filtered_features_list[[signature]]
  expr_matrix <- as.matrix(cells_with_signature[[signature]])
  if (is.null(dim(expr_matrix))) {
    expr_matrix <- matrix(expr_matrix, nrow = 1, dimnames = list(genes, colnames(counts_matrix)))
  }
  expr_df <- as.data.frame(expr_matrix)
  expr_df$Gene <- genes
  expr_df
})
presence_absence_df <- do.call(rbind, combined_list)

# Ensure that the gene column is first
presence_absence_df <- presence_absence_df %>%
  select(Gene, everything())

# Create a mapping of cell barcodes to clusters
cell_cluster_mapping <- data.frame(
  CellBarcode = colnames(counts_matrix),
  Cluster = Idents(svr.seurat.obj)[colnames(counts_matrix)]
)

# Expand presence_absence_df to have one row per gene-cell combination
expanded_df <- expand.grid(
  Gene = unique(presence_absence_df$Gene),
  CellBarcode = colnames(counts_matrix)
)

# Merge presence_absence_df with expanded_df to fill in missing combinations
merged_df <- merge(expanded_df, presence_absence_df, by = c("Gene", "CellBarcode"), all.x = TRUE)

# Fill in Presence column with FALSE where it's NA
merged_df$Presence[is.na(merged_df$Presence)] <- FALSE

# Merge with cell_cluster_mapping to get the clusters
final_df <- merge(merged_df, cell_cluster_mapping, by = "CellBarcode", all.x = TRUE)

# Write the final data frame to a file
write.table(final_df, "gene_signature_presence_absence_with_clusters.txt", sep = "\t", quote = FALSE, row.names = FALSE)








# Required libraries
library(Seurat)
library(dplyr)
library(openxlsx)

# Step 1: Read and Prepare Gene Signatures
gene_signatures <- read.table("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/DDR1hepg2_FDR025_sig_20231018_FINAL.txt", header = TRUE, sep = "\t", check.names = FALSE)

# Convert gene signatures to a list of gene sets
gene_list <- lapply(gene_signatures, as.character)
features_list <- lapply(colnames(gene_signatures), function(x) gene_signatures[[x]])
names(features_list) <- colnames(gene_signatures)

# Step 2: Get the Count Matrix
counts_matrix <- GetAssayData(svr.seurat.obj, assay = "RNA", layer = "counts")

# Step 3: Filter features_list to only include genes that are present in the counts_matrix
filtered_features_list <- lapply(features_list, function(genes) {
  intersect(genes, rownames(counts_matrix))
})

# Flatten the filtered_features_list to get a unique list of all genes
all_genes <- unique(unlist(filtered_features_list))

# Step 4: Calculate the maximum expression of each gene in each cluster using the original expression values
cluster_max_expr <- data.frame(matrix(0, nrow = length(all_genes), ncol = length(unique(Idents(svr.seurat.obj)))))
rownames(cluster_max_expr) <- all_genes
colnames(cluster_max_expr) <- levels(Idents(svr.seurat.obj))

# Populate the cluster_max_expr data frame with maximum expression values
for (gene in all_genes) {
  expr_values <- FetchData(svr.seurat.obj, vars = gene)
  expr_by_cluster <- aggregate(expr_values[[gene]], by = list(as.character(Idents(svr.seurat.obj))), FUN = max)
  cluster_max_expr[gene, expr_by_cluster$Group.1] <- expr_by_cluster$x
}

# Step 5: Determine the cluster in which each gene has the highest expression value
max_expr_df <- apply(cluster_max_expr, 1, function(row) {
  max_cluster <- which.max(row)
  max_value <- max(row)
  return(c(Cluster = colnames(cluster_max_expr)[max_cluster], Expression = max_value))
})

max_expr_df <- t(max_expr_df)
max_expr_df <- as.data.frame(max_expr_df)
max_expr_df$Gene <- rownames(max_expr_df)

max_expr_df <- max_expr_df[, c("Cluster", "Expression","Gene")]

# Step 6: Write the final data frame to an Excel file
write.xlsx(max_expr_df, "gene_signature_highest_expression_clusters.xlsx", row.names = FALSE)



##########clusters with cells
# Extract cluster information
clusters <- Idents(svr.seurat.obj)
# Create a data frame with cluster numbers and corresponding cell names
cluster_data <- data.frame(Cell = names(clusters), Cluster = as.character(clusters))
# Write the data frame to an Excel file
write.xlsx(cluster_data, "cells_in_clusters.xlsx", row.names = FALSE)



# Step 3: Filter the genes based on the gene signature
counts_matrix <- GetAssayData(svr.seurat.obj, assay = "RNA", layer =  "counts")
filtered_features_list <- lapply(features_list, function(genes) {
  intersect(genes, rownames(counts_matrix))
})
all_genes <- unique(unlist(filtered_features_list))

# Step 4: Create the expression matrix for the filtered genes
expr_matrix <- counts_matrix[all_genes, , drop = FALSE]

# Step 5: Get the cell barcodes and their clusters
clusters <- Idents(svr.seurat.obj)
cluster_info <- data.frame(Cell = names(clusters), Cluster = as.character(clusters))

# Step 6: Reorder the cell barcodes to match the order in the counts matrix
cluster_info <- cluster_info[match(colnames(counts_matrix), cluster_info$Cell), ]

# Step 7: Prepare the final output table
cluster_row <- as.character(cluster_info$Cluster)
names(cluster_row) <- colnames(expr_matrix)
final_table <- rbind(cluster_row, as.data.frame(expr_matrix))
rownames(final_table)[1] <- "Clusters"

# Add the gene names as the first column
final_table <- cbind(Gene = c("Clusters", all_genes), final_table)

# Step 8: Write the final data frame to an Excel file
write.xlsx(final_table, "gene_expression_with_clusters.xlsx", row.names = FALSE)