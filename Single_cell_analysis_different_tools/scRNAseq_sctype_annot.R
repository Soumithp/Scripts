# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# Install reticulate if not already installed
install.packages("reticulate")
library("reticulate")
reticulate::py_install("leidenalg")


setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki")

# Step 1: Read and preprocess the Seurat object (as you have done earlier)
svr_data <- Read10X_h5(filename = "C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/hg38/hg38_filtered_feature_bc_matrix.h5")
svr.seurat.obj <- CreateSeuratObject(counts = svr_data, project = "Svr", min.cells = 3, min.features = 200)
svr.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(svr.seurat.obj, pattern = "^MT-")
VlnPlot(svr.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

svr.seurat.obj <- subset(svr.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
svr.seurat.obj <- NormalizeData(svr.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
svr.seurat.obj <- FindVariableFeatures(svr.seurat.obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(svr.seurat.obj)
svr.seurat.obj <- ScaleData(svr.seurat.obj, features = all.genes)
svr.seurat.obj <- RunPCA(svr.seurat.obj, features = VariableFeatures(object = svr.seurat.obj))
# Check number of PC components
ElbowPlot(svr.seurat.obj)
svr.seurat.obj <- FindNeighbors(svr.seurat.obj, dims = 1:10)
svr.seurat.obj <- FindClusters(svr.seurat.obj, resolution = 0.8, algorithm = 3)
svr.seurat.obj <- RunUMAP(svr.seurat.obj, dims = 1:10)
DimPlot(svr.seurat.obj, reduction = "umap")





# Step 2: Read Gene Signatures
#gene_signatures <- read.table("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiver_MP_signature.txt", header = TRUE, sep = "\t", check.names = FALSE)


gene_signatures <- "C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiver_MP_allsig_wodup_sctype_format.xlsx"
tissue<- "Liver"
# prepare gene sets
gs_list <- gene_sets_prepare(gene_signatures, tissue)


# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(svr.seurat.obj[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))


# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(svr.seurat.obj[["RNA"]]$scale.data) else as.matrix(svr.seurat.obj[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. pbmc[["RNA"]]$counts or pbmc[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(svr.seurat.obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(svr.seurat.obj@meta.data[svr.seurat.obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(svr.seurat.obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# Label the clusters with their top-scoring cell type
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster == j, ]; 
  svr.seurat.obj@meta.data$sctype_classification[svr.seurat.obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Visualize the UMAP with all clusters labeled by their highest-scoring cell type
DimPlot(svr.seurat.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')


# Create a data frame containing the cell annotations
cell_annotations <- data.frame(
  Cell_Barcode = rownames(svr.seurat.obj@meta.data),
  Cluster = svr.seurat.obj@meta.data$seurat_clusters,
  Cell_Type = svr.seurat.obj@meta.data$sctype_classification
)

# Write the data frame to a tab-separated text file
write.table(cell_annotations, "C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/cell_annotations_scType_scMPliver_allsig_wodup.txt", sep = "\t", row.names = FALSE, quote = FALSE)
