

library(Seurat)
library(tidyverse)
library(hdf5r)
library(ggplot2)
library(dplyr)
library(pheatmap)
setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/")

counts <- Read10X_h5("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/hg38/hg38_filtered_feature_bc_matrix.h5")
df <- CreateSeuratObject(counts = counts)
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
df <- subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
df <- NormalizeData(df)
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(df)
df <- ScaleData(df, features = all.genes)
df <- RunPCA(df, features = VariableFeatures(object = df))
# df <- FindNeighbors(df)
# df <- FindClusters(df, resolution = 1.8)
# df <- RunUMAP(df, dims = 1:10)
# DimPlot(df, reduction = "umap", pt.size = 3)



#Using PCA embeddings for clustering (more common for dimensionality reduction)
pca_data <- as.matrix(df@reductions$pca@cell.embeddings)
# Calculate the distance matrix on the chosen data
distance_matrix <- dist(pca_data)  # Use PCA embeddings for hierarchical clustering
# Perform hierarchical clustering
hclust_result <- hclust(distance_matrix, method = "ward.D2")
# Cut the tree to form clusters
clusters <- cutree(hclust_result, k = 10)  # Adjust 'k' to your desired number of clusters
# Add the clusters to the Seurat object metadata
df$hclust_clusters <- as.factor(clusters)
# Extract the UMAP coordinates
umap_data <- as.data.frame(df@reductions$umap@cell.embeddings)
# Plot the UMAP with hierarchical clusters
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = df$hclust_clusters)) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "UMAP Plot with Hierarchical Clusters", color = "Cluster") +
  theme_minimal()

# Use the PCA data for the heatmap
pca_data <- as.matrix(df@reductions$pca@cell.embeddings)
# Order the data based on hierarchical clustering
ordered_pca_data <- pca_data[order(df$hclust_clusters), ]
# Create a heatmap
pheatmap(
  ordered_pca_data, 
  scale = "row",
  cluster_rows = TRUE,   # Rows are already ordered by clustering
  cluster_cols = TRUE,   # Columns can be left as is or clustered depending on your needs
  show_rownames = FALSE,  # Disable row names for a cleaner heatmap
  show_colnames = TRUE,  # Disable column names for a cleaner heatmap
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap of PCA Embeddings by Hierarchical Clusters", 
  cellheight = 0.5, # Adjust this value to make the rows more compact
  cellwidth = 10
)


scLiver_MP_signature <- read.delim("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiverMP_sig_allcombined_wnodup.txt")

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

view(df)

DimPlot(df, reduction = "umap", group.by = "cell_type", pt.size = 3)

# Extract the module scores from the Seurat object metadata
module_scores <- df@meta.data[, grep("^Celltype_", colnames(df@meta.data))]

# Transpose the data for heatmap plotting
module_scores_t <- t(module_scores)

# Define a custom color palette: low values in blue, high values in red
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Generate the heatmap
pheatmap(
  module_scores_t, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  scale = "row", 
  show_rownames = TRUE, 
  show_colnames = FALSE, 
  main = "Heatmap of Module Scores",
  color = color_palette,
  cellheight = 20, # Adjust this value to make the rows more compact
  cellwidth = 1
)

####heatmap with orginal values
pheatmap(
  module_scores_t, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  show_rownames = TRUE, 
  show_colnames = FALSE, 
  main = "Heatmap of Module Scores",
  color = color_palette,
  cellheight = 20, # Adjust this value to make the rows more compact
  cellwidth = 1
)


# Order module scores by hierarchical clusters
ordered_module_scores <- module_scores[order(df$hclust_clusters), ]
# Transpose the data for heatmap plotting
ordered_module_scores_t <- t(ordered_module_scores)
# Generate the heatmap for ordered module scores
pheatmap(
  ordered_module_scores_t, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  scale = "row", 
  show_rownames = TRUE, 
  show_colnames = FALSE, 
  main = "Heatmap of Module Scores Ordered by Hierarchical Clustering",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cellheight = 20, # Adjust this value to make the rows more compact
  cellwidth = 1
)

# ####for clusters and score of individual cells
# # Extract the cluster information as numeric values
# cluster_info <- as.numeric(df$seurat_clusters)  # Convert to numeric if needed
# # Combine the cluster information with the module scores
# heatmap_data <- cbind(Cluster = cluster_info, module_scores)
# # Order cells by their cluster
# heatmap_data <- heatmap_data[order(heatmap_data$Cluster), ]
# # Exclude the cluster information before transposing
# module_scores_only <- heatmap_data[, -1]  # Remove the Cluster column
# # Transpose the module scores for heatmap plotting
# module_scores_t <- t(module_scores_only)  # Transpose only the numeric data
# # Set the row names to the numerical cluster values
# row.names(module_scores_t) <- cluster_info
# # Generate the heatmap
# pheatmap(
#   module_scores_t, 
#   cluster_rows = TRUE, 
#   cluster_cols = FALSE,  # Disable clustering of cells to maintain cluster order
#   scale = "row", 
#   show_rownames = TRUE, 
#   show_colnames = TRUE,   # Show columns since each column is a cell
#   main = "Heatmap of Module Scores by Cluster"
# )

# Specify the cell type to include
selected_cell_type <- "Hepatocyte"
# Subset the Seurat object to include only the specified cell type
subset_df <- subset(df, subset = cell_type == selected_cell_type)
# Extract the specified columns from the metadata
selected_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt")
# Create a new metadata data frame with only the selected columns
new_metadata <- subset_df@meta.data[, selected_columns]
# Create a new Seurat object using the subsetted data and the new metadata
subset_df_new <- CreateSeuratObject(
  counts = GetAssayData(subset_df, slot = "counts"),
  meta.data = new_metadata
)
# Add the UMAP embeddings to the new Seurat object
subset_df_new[["umap"]] <- subset_df[["umap"]]
# Check the metadata of the new Seurat object
head(subset_df_new@meta.data)

subset_df_new <- NormalizeData(subset_df_new)
all.genes <- rownames(subset_df_new)
subset_df_new <- FindVariableFeatures(subset_df_new, selection.method = "vst", nfeatures = 1000)
subset_df_new <- ScaleData(subset_df_new, features = all.genes)
subset_df_new <- RunPCA(subset_df_new, features = VariableFeatures(object = subset_df_new))
subset_df_new <- FindNeighbors(subset_df_new)
subset_df_new <- FindClusters(subset_df_new, resolution = 1.5)
subset_df_new <- RunUMAP(subset_df_new, dims = 1:10)
DimPlot(subset_df_new, reduction = "umap", pt.size = 3)


scLiverMP_hepa_sig <- read.delim("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiverMP_hepatocyte_sig.txt")


###hierarchial clustering
#Using PCA embeddings for clustering (more common for dimensionality reduction)
pca_data <- as.matrix(subset_df_new@reductions$pca@cell.embeddings)
# Calculate the distance matrix on the chosen data
distance_matrix <- dist(pca_data)  # Use PCA embeddings for hierarchical clustering
# Perform hierarchical clustering
hclust_result <- hclust(distance_matrix, method = "ward.D2")
# Cut the tree to form clusters
clusters <- cutree(hclust_result, k = 6)  # Adjust 'k' to your desired number of clusters
# Add the clusters to the Seurat object metadata
subset_df_new$hclust_clusters <- as.factor(clusters)
# Extract the UMAP coordinates
umap_data <- as.data.frame(subset_df_new@reductions$umap@cell.embeddings)
# Plot the UMAP with hierarchical clusters
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = subset_df_new$hclust_clusters)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "pink")) +
  labs(title = "UMAP Plot with Hierarchical Clusters", color = "Cluster") +
  theme_minimal()



########using heatmap_for hepatocyte annotation
# Use the PCA data for the heatmap
pca_data <- as.matrix(subset_df_new@reductions$pca@cell.embeddings)
# Order the data based on hierarchical clustering
ordered_pca_data <- pca_data[order(subset_df_new$hclust_clusters), ]
# Create a heatmap
pheatmap(
  ordered_pca_data, 
  scale = "row",
  cluster_rows = TRUE,   # Rows are already ordered by clustering
  cluster_cols = TRUE,   # Columns can be left as is or clustered depending on your needs
  show_rownames = FALSE,  # Disable row names for a cleaner heatmap
  show_colnames = TRUE,  # Disable column names for a cleaner heatmap
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap of PCA Embeddings by Hierarchical Clusters", 
  cellheight = 1, # Adjust this value to make the rows more compact
  cellwidth = 10
)










signatures_hepatocytes = list()

for(cell_type in colnames(scLiverMP_hepa_sig)){
  print(paste('Extracting the cell type:',cell_type,'...'))
  signature_list = scLiverMP_hepa_sig[[cell_type]]
  signatures_hepatocytes[[cell_type]] = signature_list[which(signature_list!="")]
}

for(i in 1:ncol(scLiverMP_hepa_sig)){
  print(i)
  subset_df_new <- AddModuleScore(object = subset_df_new, features = list(signatures_hepatocytes[[i]]), name = paste0('Celltype_',names(signatures_hepatocytes)[i]))
}

scores= subset_df_new@meta.data[, grep("^Celltype_", colnames(subset_df_new@meta.data))]
head(scores)
subset_df_new$cell_type = apply(scores, 1, function(x) {
  cell_types <- names(signatures_hepatocytes)
  cell_types[which.max(x)]
})



DimPlot(subset_df_new, reduction = "umap", group.by = "cell_type", pt.size = 3)

