setwd("C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki")
# Step 1: Read and preprocess the Seurat object
svr_data <- Read10X_h5(filename = "C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/hg38/hg38_filtered_feature_bc_matrix.h5")
svr.seurat.obj <- CreateSeuratObject(counts = svr_data, project = "Svr", min.cells = 3, min.features = 200)
svr.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(svr.seurat.obj, pattern = "^MT-")
svr.seurat.obj <- subset(svr.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
svr.seurat.obj <- NormalizeData(svr.seurat.obj)
svr.seurat.obj <- FindVariableFeatures(svr.seurat.obj, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(svr.seurat.obj)
svr.seurat.obj <- ScaleData(svr.seurat.obj, features = all.genes)
svr.seurat.obj <- RunPCA(svr.seurat.obj, features = VariableFeatures(object = svr.seurat.obj))
svr.seurat.obj <- FindNeighbors(svr.seurat.obj, dims = 1:10)
svr.seurat.obj <- FindClusters(svr.seurat.obj, resolution = 0.7)
svr.seurat.obj <- RunUMAP(svr.seurat.obj, dims = 1:10)
DimPlot(svr.seurat.obj, reduction = "umap", pt.size = 4)



# Load necessary libraries
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(celldex)  # Contains pre-built references
library(ggplot2)
# Step 1: Preprocess the Seurat object (as you already did)
# Assuming svr.seurat.obj is already created and preprocessed

# Step 2: Convert Seurat object to SingleCellExperiment object for SingleR
sce <- as.SingleCellExperiment(svr.seurat.obj)

# Step 3: Load a reference dataset for liver cells or use pre-built references from celldex
# Option 1: Using pre-built reference data (e.g., Human Primary Cell Atlas)
ref <- celldex::HumanPrimaryCellAtlasData()

# Option 2: Use a liver-specific reference dataset if you have one:
# ref <- readRDS("path_to_liver_reference_dataset.rds")

# Step 4: Run SingleR for cluster-based annotation
# Note: cluster-based annotation helps assign a cell type label to each cluster
clusters <- svr.seurat.obj@meta.data$seurat_clusters

# Run SingleR
singleR_results <- SingleR(test = sce, ref = ref, labels = ref$label.main, clusters = clusters)

# Step 5: Add the SingleR annotations back to the Seurat object
svr.seurat.obj$SingleR.labels <- singleR_results$labels[match(svr.seurat.obj$seurat_clusters, rownames(singleR_results))]

# Step 6: Visualize the results
# You can visualize the UMAP plot and color by SingleR-annotated cell types
DimPlot(svr.seurat.obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, pt.size = 4) +
  ggtitle("UMAP with SingleR Cell Type Annotations")

# Optionally, check the distribution of cell types across clusters
table(svr.seurat.obj$SingleR.labels, svr.seurat.obj$seurat_clusters)

# Step 7: Save the annotated object for further use
saveRDS(svr.seurat.obj, file = "svr_seurat_with_singleR_annotations.rds")



######genesets extraction
library(celldex)
ref <- celldex::HumanPrimaryCellAtlasData()
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
