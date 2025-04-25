
library(Seurat)
library(tidyverse)
library(hdf5r)
library(ggplot2)
library(dplyr)
library(pheatmap)
setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/")

scLiver_MP_signature <- read.delim("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiverMP_sig_allcombined_wnodup.txt")
counts <- Read10X_h5("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/hg38/hg38_filtered_feature_bc_matrix.h5")
df <- CreateSeuratObject(counts = counts)
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
df <- subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
df <- NormalizeData(df)
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(df)
df <- ScaleData(df, features = all.genes)
df <- RunPCA(df, features = VariableFeatures(object = df))
df <- FindNeighbors(df)
df <- FindClusters(df, resolution = 1.8)
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

view(df)

DimPlot(df, reduction = "umap", group.by = "cell_type", pt.size = 3)


cluster1 <- df@meta.data %>%
  filter(seurat_clusters == 7)




# ###########looking into hepatocytes
# #Specify the clusters to include
# selected_clusters <- c(0, 1, 2, 3, 7)
# # Subset the Seurat object to include only the specified clusters
# subset_df <- subset(df, idents = selected_clusters)

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


scLiverMP_hepa_sig <- read.delim("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/DDR1hepg2_FDR025_sig_20231018_FINAL.txt")


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

###difference betwen up and down
scores$Difference <- scores$Celltype_DDR1hepg2_FDR025_UP1 - scores$Celltype_DDR1hepg2_FDR025_DN1
subset_df_new$Difference <- scores$Difference


# Add the `up`, `down`, and `difference` columns to the Seurat object metadata
subset_df_new$Up <- scores$Celltype_DDR1hepg2_FDR025_UP1
subset_df_new$Down <- scores$Celltype_DDR1hepg2_FDR025_DN1
subset_df_new$Difference <- scores$Difference

# Create individual UMAP plots with red color points
p_up <- FeaturePlot(subset_df_new, features = "Up", pt.size = 4, cols = c("lightgrey", "red")) +
  ggplot2::ggtitle("UMAP plot for Up") +
  NoLegend()

p_down <- FeaturePlot(subset_df_new, features = "Down", pt.size = 4, cols = c("lightgrey", "red")) +
  ggplot2::ggtitle("UMAP plot for Down") +
  NoLegend()

p_difference <- FeaturePlot(subset_df_new, features = "Difference", pt.size = 4, cols = c("lightgrey", "red")) +
  ggplot2::ggtitle("UMAP plot for Difference") +
  NoLegend()

# Combine the UMAP plots into one frame with adjusted plot dimensions
combined_plot <- plot_grid(
  p_up, p_down, p_difference, 
  ncol = 3, 
  align = 'hv',
  rel_heights = c(0.2),  # Adjust the relative height of the plots
  rel_widths = c(1, 1, 1)  # Adjust the relative width of the plots
)

# Print the combined plot
print(combined_plot)




