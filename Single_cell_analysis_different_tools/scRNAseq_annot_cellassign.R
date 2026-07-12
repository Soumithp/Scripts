BiocManager::install("cellassign")


library(Seurat)
library(SeuratObject)
library(cellassign)
library(tensorflow)
library(openxlsx)
library(dplyr)
library(devtools)
library(scran)

install.packages("devtools") # If not already installed
devtools::install_github("Irrationone/cellassign")

install.packages("reticulate")

install.packages("tensorflow")
library(tensorflow)
install_tensorflow(extra_packages = "tensorflow-probability")
tensorflow::tf_config()

tensorflow::install_tensorflow()



BiocManager::install('cellassign')


install.packages("glue")

setwd("C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki")

source("cellassign.R")
source("utils.R")
source("simulate.R")
source("inference-tensorflow.R")

# Set working directory
setwd("C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki")

# Step 1: Read and preprocess the Seurat object
svr_data <- Read10X_h5(filename = "C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/hg38/hg38_filtered_feature_bc_matrix.h5")
svr.seurat.obj <- CreateSeuratObject(counts = svr_data, project = "Svr", min.cells = 3, min.features = 200)
svr.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(svr.seurat.obj, pattern = "^MT-")
VlnPlot(svr.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

svr.seurat.obj <- subset(svr.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
svr.seurat.obj <- NormalizeData(svr.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
svr.seurat.obj <- FindVariableFeatures(svr.seurat.obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(svr.seurat.obj)
svr.seurat.obj <- ScaleData(svr.seurat.obj, features = all.genes)
svr.seurat.obj <- RunPCA(svr.seurat.obj, features = VariableFeatures(object = svr.seurat.obj))
ElbowPlot(svr.seurat.obj)
svr.seurat.obj <- FindNeighbors(svr.seurat.obj, dims = 1:10)
svr.seurat.obj <- FindClusters(svr.seurat.obj, resolution = 0.8, algorithm = 3)
svr.seurat.obj <- RunUMAP(svr.seurat.obj, dims = 1:10)
DimPlot(svr.seurat.obj, reduction = "umap")



# Step 1: Read the gene signatures file
gene_signatures <- read.table("C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiverMP_sig_allcombined_wnodup.txt", 
                              header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# Step 2: Prepare a binary marker gene matrix
all_marker_genes <- unique(na.omit(unlist(gene_signatures)))
all_marker_genes <- all_marker_genes[all_marker_genes != ""]  # Remove empty strings

# Initialize the marker gene matrix
marker_gene_matrix <- matrix(0, nrow = length(all_marker_genes), ncol = ncol(gene_signatures))
rownames(marker_gene_matrix) <- all_marker_genes
colnames(marker_gene_matrix) <- colnames(gene_signatures)

# Populate the marker gene matrix
for (cell_type in colnames(gene_signatures)) {
  genes <- na.omit(gene_signatures[[cell_type]])
  genes <- genes[genes != ""]
  marker_gene_matrix[genes, cell_type] <- 1
}

# Step 3: Get common genes between marker genes and the expression matrix (sce)
common_genes <- intersect(rownames(marker_gene_matrix), rownames(svr.seurat.obj))

# Filter the marker gene matrix
marker_gene_matrix_filtered <- marker_gene_matrix[common_genes, , drop = FALSE]

# Filter the Seurat object expression data
sce_filtered <- svr.seurat.obj[common_genes, , drop = FALSE]

# Step 4: Convert the filtered Seurat object to a SingleCellExperiment object
# Ensure you're extracting the "counts" data from Seurat and converting it into an SCE object
counts_matrix <- GetAssayData(svr.seurat.obj, slot = "counts")[common_genes, ]
sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

# Step 5: Compute Size Factors using scran
# This will normalize the data and provide size factors for each cell
sce <- computeSumFactors(sce)

# Step 6: Run CellAssign
fit <- cellassign(
  exprs_obj = sce,  # SingleCellExperiment object
  marker_gene_info = marker_gene_matrix_filtered,  # Filtered marker gene matrix
  s = sizeFactors(sce),  # Use the size factors computed by scran
  learning_rate = 1e-2,   # Adjust learning rate if needed
  shrinkage = TRUE        # Enable shrinkage for better accuracy
)

# Step 7: Add CellAssign results to Seurat object and visualize
svr.seurat.obj$cellassign_labels <- fit$cell_type

# Visualize UMAP with CellAssign labels
DimPlot(svr.seurat.obj, reduction = "umap", group.by = "cellassign_labels", label = TRUE, repel = TRUE)

# Export annotated cells
cell_annotations <- data.frame(
  Cell_Barcode = rownames(svr.seurat.obj@meta.data),
  Cluster = svr.seurat.obj@meta.data$seurat_clusters,
  Cell_Type = svr.seurat.obj$cellassign_labels
)

write.table(cell_annotations, "C:/Users/soumith/OneDrive - University of Texas Southwestern/Desktop/cell_annotations_cellassign.txt", sep = "\t", row.names = FALSE, quote = FALSE)