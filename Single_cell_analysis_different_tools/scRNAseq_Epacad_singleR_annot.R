setwd("C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/scRNAseq/organised_2024.09.30")
# Step 1: Read and preprocess the Seurat object
svr_data <- Read10X_h5(filename = "C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/scRNAseq/organised_2024.09.30/hg38_filtered_feature_bc_matrix.h5")
svr.seurat.obj <- CreateSeuratObject(counts = svr_data, project = "Svr", min.cells = 3, min.features = 200)
svr.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(svr.seurat.obj, pattern = "^MT-")
svr.seurat.obj <- subset(svr.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
svr.seurat.obj <- NormalizeData(svr.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
svr.seurat.obj <- FindVariableFeatures(svr.seurat.obj, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(svr.seurat.obj)
svr.seurat.obj <- ScaleData(svr.seurat.obj, features = all.genes)
svr.seurat.obj <- RunPCA(svr.seurat.obj, features = VariableFeatures(object = svr.seurat.obj))
svr.seurat.obj <- FindNeighbors(svr.seurat.obj, dims = 1:10)
svr.seurat.obj <- FindClusters(svr.seurat.obj, resolution = 0.8, algorithm = 3)
svr.seurat.obj <- RunUMAP(svr.seurat.obj, dims = 1:10)
DimPlot(svr.seurat.obj, reduction = "umap", pt.size = 2)



# Load necessary libraries
library(reticulate)
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(celldex)  # Contains pre-built references
library(ggplot2)
############singleR usign each cell type annotations using celldex::MonacoImmuneData()
# sce <- as.SingleCellExperiment(svr.seurat.obj)
# ref <- celldex::MonacoImmuneData()
# singleR_results <- SingleR(test = sce, ref = ref, labels = ref$label.main)
# svr.seurat.obj$SingleR.labels <- singleR_results$labels
# DimPlot(svr.seurat.obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, pt.size = 2) +
#   ggtitle("UMAP with SingleR Cell Type Annotations")
# 
# table(svr.seurat.obj$SingleR.labels, svr.seurat.obj$seurat_clusters)
# saveRDS(svr.seurat.obj, file = "svr_seurat_with_singleR_annotations.rds")

############singleR usign each cell type annotations using humanprimarycellatlasdata
sce <- as.SingleCellExperiment(svr.seurat.obj)
ref <-  celldex::HumanPrimaryCellAtlasData()
singleR_results <- SingleR(test = sce, ref = ref, labels = ref$label.main)
svr.seurat.obj$SingleR.labels <- singleR_results$labels
DimPlot(svr.seurat.obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, pt.size = 2) +
  ggtitle("UMAP with SingleR Cell Type Annotations")

table(svr.seurat.obj$SingleR.labels, svr.seurat.obj$seurat_clusters)
saveRDS(svr.seurat.obj, file = "svr_seurat_with_singleR_annotations.rds")

# Save the cell type-based annotations
cell_type_based_annotations <- data.frame(
  Cell_Barcode = rownames(svr.seurat.obj@meta.data),  # Barcodes of cells
  Cluster = svr.seurat.obj$seurat_clusters,  # Cluster IDs for each cell
  SingleR_Label = svr.seurat.obj$SingleR.labels  # SingleR labels for each cell
)

# Write the cell type-based annotations to a text file
write.table(cell_type_based_annotations, "cell_type_based_annotations.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ########using singleR clustering based annotations using monacoimmunedata
# sce <- as.SingleCellExperiment(svr.seurat.obj)
# ref <- celldex::MonacoImmuneData()
# clusters <- svr.seurat.obj@meta.data$seurat_clusters
# singleR_results <- SingleR(test = sce, ref = ref, labels = ref$label.main, clusters = clusters)
# svr.seurat.obj$SingleR.labels <- singleR_results$labels[match(svr.seurat.obj$seurat_clusters, rownames(singleR_results))]
# DimPlot(svr.seurat.obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, pt.size = 2) +
#   ggtitle("UMAP with SingleR Cluster-Based Annotations")
# 
# table(svr.seurat.obj$SingleR.labels, svr.seurat.obj$seurat_clusters)
# saveRDS(svr.seurat.obj, file = "svr_seurat_with_singleR_cluster_annotations.rds")
# 
# ########using singleR clustering based annotations using celldex::BlueprintEncodeData()
# sce <- as.SingleCellExperiment(svr.seurat.obj)
# ref <- celldex::BlueprintEncodeData()
# clusters <- svr.seurat.obj@meta.data$seurat_clusters
# singleR_results <- SingleR(test = sce, ref = ref, labels = ref$label.main, clusters = clusters)
# svr.seurat.obj$SingleR.labels <- singleR_results$labels[match(svr.seurat.obj$seurat_clusters, rownames(singleR_results))]
# DimPlot(svr.seurat.obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, pt.size = 2) +
#   ggtitle("UMAP with SingleR Cluster-Based Annotations")
# 
# table(svr.seurat.obj$SingleR.labels, svr.seurat.obj$seurat_clusters)
# saveRDS(svr.seurat.obj, file = "svr_seurat_with_singleR_cluster_annotations.rds")

########using singleR clustering based annotations using celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(svr.seurat.obj)
ref <- celldex::HumanPrimaryCellAtlasData()
clusters <- svr.seurat.obj@meta.data$seurat_clusters
singleR_results <- SingleR(test = sce, ref = ref, labels = ref$label.main, clusters = clusters)
svr.seurat.obj$SingleR.labels <- singleR_results$labels[match(svr.seurat.obj$seurat_clusters, rownames(singleR_results))]
DimPlot(svr.seurat.obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, pt.size = 2) +
  ggtitle("UMAP with SingleR Cluster-Based Annotations")

table(svr.seurat.obj$SingleR.labels, svr.seurat.obj$seurat_clusters)
saveRDS(svr.seurat.obj, file = "svr_seurat_with_singleR_cluster_annotations.rds")


############annotations
# Save the cluster-based annotations
cluster_based_annotations <- data.frame(
  Cluster = clusters,  # The cluster IDs
  Cell_Barcode = rownames(svr.seurat.obj@meta.data),  # Barcodes of cells
  SingleR_Label = svr.seurat.obj$SingleR.labels  # Cluster-based SingleR labels
)
# Write the annotations to a text file
write.table(cluster_based_annotations, "SingleR_cluster_based_annotations.txt", sep = "\t", row.names = FALSE, quote = FALSE)


























######genesets extraction
library(celldex)
ref <- celldex::BlueprintEncodeData()
# View the available metadata for the reference
colData(ref)

# View the matrix of gene expression profiles for each cell type
expression_matrix <- assay(ref, "logcounts")

# View the cell type labels in the reference
head(ref$label.main)

# Filter for Hepatocytes and T Cells in the reference
hepatocyte_idx <- which(cell_types == "Hepatocytes")
tcell_idx <- which(grepl("T_cell", cell_types, ignore.case = TRUE))

# Extract gene expression for Hepatocytes and T Cells
hepatocyte_expression <- expression_matrix[, hepatocyte_idx]
tcell_expression <- expression_matrix[, tcell_idx]

# Calculate the average expression for Hepatocytes and T Cells
avg_hepatocyte_expr <- rowMeans(hepatocyte_expression)
avg_tcell_expr <- rowMeans(tcell_expression)

# Sort the genes by expression level
all_hepatocyte_genes <- names(sort(avg_hepatocyte_expr, decreasing = TRUE))
all_tcell_genes <- names(sort(avg_tcell_expr, decreasing = TRUE))

# Output all genes for Hepatocytes and T Cells
print("All Hepatocyte genes:")
print(all_hepatocyte_genes)

print("All T Cell genes:")
print(all_tcell_genes)

# Save the gene lists to CSV files
write.csv(all_hepatocyte_genes, "all_hepatocyte_genes.csv", row.names = FALSE)
write.csv(all_tcell_genes, "all_tcell_genes.csv", row.names = FALSE)
