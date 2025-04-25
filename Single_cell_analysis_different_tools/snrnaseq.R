# RStudio 2023.09.0+463 "Desert Sunflower" Release (b51c81cc303d4b52b010767e5b30438beb904641, 2023-09-25) for windows
# Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) RStudio/2023.09.0+463 Chrome/114.0.5735.289 Electron/25.5.0 Safari/537.36
setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/scRNAseq")

library(Seurat)
library(tidyverse)
library(hdf5r)
library(ggplot2)

#for annotation
library(SingleR)
library(celldex)

#to make the excel format matrix
veh_data_df <- as.data.frame(as.matrix(Aggr_data))
write.csv(veh_data_df, "snRNAseq_Aggre.csv", row.names = TRUE)

#reading the filtered data
veh_data<- Read10X_h5(filename= "C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/scRNAseq/Veh1/outs/filtered_feature_bc_matrix.h5")

#creating suerat object(non-normalized data)
veh.seurat.obj<- CreateSeuratObject(counts= veh_data,  project = "Veh", min.cells= 3, min.features= 200)
# str(veh.seurat.obj)

##Qc %mt reads
veh.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(veh.seurat.obj, pattern = "Mt")
VlnPlot(veh.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(veh.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering
veh.seurat.obj <- subset(veh.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 12500 & percent.mt < 5)


##normalizing data
# 3. Normalize data 
veh.seurat.obj <- NormalizeData(veh.seurat.obj)




library(Seurat)
veh.seurat.obj <- subset(veh.seurat.obj, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)
# Step 3: Normalize data
veh.seurat.obj <- NormalizeData(veh.seurat.obj)
# Step 4: Extract the normalized data matrix
normalized_matrix <- GetAssayData(veh.seurat.obj, layer = "data")
# Step 5: Save the matrix to a text file
write.table(as.matrix(normalized_matrix), file = "veh1_filtered_normalized_matrix.txt", sep = "\t", quote = FALSE, col.names = NA)
print("Matrix saved to 'filtered_normalized_matrix.txt'")


# ###saving normalized data
# # Extract normalized data
# normalized_data <- GetAssayData(veh.seurat.obj, slot = "data")
# normalized_data <- as.data.frame(as.matrix(normalized_data))
# # Write to CSV
# write.csv(normalized_data, "snRNAseq_Aggr_norm.csv", row.names = TRUE)
#to see the commands which were used before normalizaton 
#@commands in the str(veh.seurat.obj)

# 4.Identify highly variable features 
veh.seurat.obj <- FindVariableFeatures(veh.seurat.obj, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(veh.seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(veh.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scaling
all.genes <- rownames(veh.seurat.obj)
veh.seurat.obj <- ScaleData(veh.seurat.obj, features = all.genes)

# 6. Perform Linear dimensionality reduction
veh.seurat.obj <- RunPCA(veh.seurat.obj, features = VariableFeatures(object = veh.seurat.obj))

# visualize PCA results
print(veh.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(veh.seurat.obj, dims = 1:5, cells = 500, balanced = TRUE)
# determine dimensionality of the data
ElbowPlot(veh.seurat.obj)


# 7. Clustering
veh.seurat.obj <- FindNeighbors(veh.seurat.obj, dims = 1:17)

# understanding resolution
veh.seurat.obj <- FindClusters(veh.seurat.obj, resolution = 1.4)

DimPlot(veh.seurat.obj, group.by = "RNA_snn_res.1.4", label = TRUE)

# setting identity of clusters
Idents(veh.seurat.obj) <- "RNA_snn_res.1.4"
Idents(veh.seurat.obj)


veh.seurat.obj <- RunUMAP(veh.seurat.obj, dims = 1:17)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(veh.seurat.obj, reduction = "umap")


veh.seurat.obj@meta.data


##############################################addmodulescore##########
rownames(veh.seurat.obj) <- toupper(rownames(veh.seurat.obj))
head(rownames(veh.seurat.obj))

# Assuming you have read your gene signatures into `gene_signatures` DataFrame
gene_signatures <- read.table("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/scRNAseq/msZone_signature.txt", header = TRUE, sep= "\t", check.names = FALSE)

# Convert the entire data frame to uppercase using base R
gene_signatures <- data.frame(lapply(gene_signatures, toupper), stringsAsFactors = FALSE)



View(gene_signatures)
# Transpose and convert to list
gene_list <- lapply(gene_signatures, as.character)

# Example: assuming gene_signatures is already a data frame with genes in columns under each cell type
features_list <- lapply(colnames(gene_signatures), function(x) gene_signatures[[x]])
names(features_list) <- colnames(gene_signatures)

features_list

# AddModuleScore for each cell type
for (cell_type in names(features_list)) {
  # Create a concise score name based on the cell type name
  score_name <- gsub("[[:punct:]]|\\s+", "_", cell_type)  # Replace punctuation and spaces with underscores
  # Apply AddModuleScore for each set of genes associated with each cell type
  # Ensure the features for the module are wrapped in a list
  veh.seurat.obj <- AddModuleScore(object = veh.seurat.obj, features = list(features_list[[cell_type]]), name = score_name)
}

# Extract the added score column names from the metadata
view(veh.seurat.obj)
print(colnames(veh.seurat.obj@meta.data))

score_columns <- grep("msZone", colnames(veh.seurat.obj@meta.data), value = TRUE)
print("Score columns identified:")
print(score_columns)

# Function to create feature plots in groups of 10
create_feature_plots <- function(seurat_obj, features, group_size = 10) {
  n <- length(features)
  print(paste("Number of features to plot:", n))
  groups <- split(features, ceiling(seq_along(features) / group_size))
  
  plot_list <- list()
  
  for (i in seq_along(groups)) {
    print(paste("Creating plot for group", i))
    plot <- FeaturePlot(seurat_obj, features = groups[[i]], ncol = 5)  # Adjust ncol if needed
    plot_list[[i]] <- plot  # Store each plot in a list
  }
  
  return(plot_list)
}

# Generate the feature plots
feature_plots <- create_feature_plots(veh.seurat.obj, score_columns, group_size = 10)

# Check if plots were generated
if (length(feature_plots) == 0) {
  stop("No plots were generated. Please check the feature columns and data.")
}

# Print the first plot to ensure it is generated correctly
print(feature_plots[[1]])

# Save all the plots to files
output_dir <- "plots"  # Define the output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (i in seq_along(feature_plots)) {
  # Define the filename
  filename <- file.path(output_dir, paste0("feature_plot_", i, ".jpeg"))
  # Save the plot using ggsave
  ggsave(filename = filename, plot = feature_plots[[i]], width = 25, height = 11 )
  # Print a message confirming the save
  cat("Saved plot to", filename, "\n")
}


# Find the average module score per cluster
average_scores <- AverageExpression (veh.seurat.obj, group.by= "seurat_clusters" , return.seurat = TRUE)

# You might need to determine which score corresponds to which cluster manually or automate by finding the maximum
for (i in names(features_list)) {
  cluster <- which.max(average_scores@meta.data[[paste0(i, '1')]])
  Idents(veh.seurat.obj, cells = WhichCells(veh.seurat.obj, idents = cluster)) <- i
}


library(dplyr)


# Extract metadata to a data frame
metadata_df <- veh.seurat.obj@meta.data

# Assume score columns follow the naming pattern 'score_'
score_columns <- grep("msZone", colnames(metadata_df), value = TRUE)

# Create a new column in this data frame that contains the name of the score column with the highest value
metadata_df$highest_score_label <- apply(metadata_df[, score_columns], 1, function(x) {
  names(x)[which.max(x)]
})

# Define Mode function since R does not have a built-in mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Aggregate this information at the cluster level
cluster_summary <- metadata_df %>%
  group_by(seurat_clusters = metadata_df$seurat_clusters) %>%
  summarise(MostCommonScore = Mode(highest_score_label)) 


# Merge this summary back to the original metadata
metadata_df$cluster_label <- cluster_summary$MostCommonScore[match(metadata_df$seurat_clusters, cluster_summary$seurat_clusters)]

# Put the updated metadata back into the Seurat object
veh.seurat.obj <- AddMetaData(veh.seurat.obj, metadata = metadata_df$cluster_label, col.name = "cluster_label")

veh.seurat.obj@meta.data

# Plot UMAP with labeled clusters based on the highest score
DimPlot(veh.seurat.obj, reduction = "umap", group.by = "cluster_label", label = TRUE)


##########tsne map

# Check if t-SNE is already computed
if (!"tsne" %in% names(veh.seurat.obj@reductions)) {
  # Run t-SNE
  veh.seurat.obj <- RunTSNE(veh.seurat.obj, dims = 1:17)  # Adjust dimensions as necessary
}

# Check to confirm t-SNE is there
print(names(veh.seurat.obj@reductions))

# Plot t-SNE with labeled clusters based on the highest score
DimPlot(veh.seurat.obj, reduction = "tsne", group.by = "cluster_label", label = TRUE, repel = TRUE)






#########extracting the cell annotations
# Assuming 'veh.seurat.obj' is your Seurat object
metadata_df <- veh.seurat.obj@meta.data

# Extract only the cell barcodes and cluster_label
cell_types_df <- data.frame(
  cell_barcode = rownames(metadata_df),
  cluster_label = metadata_df$cluster_label
)

# View the first few rows to check
head(cell_types_df)

cell_types_matrix <- as.matrix(cell_types_df)

# Set the row names as cell barcodes for the matrix
rownames(cell_types_matrix) <- cell_types_df$cell_barcode

# Optionally remove the barcode column if it's redundant in row names
cell_types_matrix <- cell_types_matrix[, "cluster_label", drop = FALSE]

# Check the matrix
print(cell_types_matrix)

# Write the data frame to a CSV file
write.csv(cell_types_df, "Aggre_msZone_annot.csv", row.names = FALSE)

# Or write the matrix if preferred
#write.csv(cell_types_matrix, "Veh1_annot.csv", row.names = TRUE)


###########annotation using sccatch
install.packages("scCATCH")
library(scCATCH)

# Load your gene signatures
signature_data <- read.table("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/scRNAseq/msZone_signature.txt", header = TRUE, check.names = FALSE, row.names = 1)

# Convert this data frame to a list where each element is a vector of gene names for a cell type
gene_signatures <- lapply(signature_data, as.character)

# Make sure that each column becomes an element in the list
gene_signatures <- lapply(colnames(signature_data), function(x) as.character(signature_data[,x]))
names(gene_signatures) <- colnames(signature_data)

# Check the structure of gene_signatures to ensure it is correct
str(gene_signatures)

custom_marker= data.frame(gene_signatures)

# Ensure row names of Seurat object are all uppercase (common practice for gene names)
rownames(veh.seurat.obj) <- toupper(rownames(veh.seurat.obj))

# Extract the count matrix from the Seurat object
count_matrix <- GetAssayData(veh.seurat.obj, slot = "counts")
# Extract cluster information from the Seurat object
cluster_info <- veh.seurat.obj$seurat_clusters

# Convert cluster_info to a factor if it's not already
cluster_info <- as.character(cluster_info)

# Create a scCatch object using the count matrix
scCatch_obj <- scCATCH::createscCATCH(data = count_matrix, cluster =cluster_info )

str(scCatch_obj)


scCatch_obj <- findmarkergene(scCatch_obj, cluster = cluster_info, if_use_custom_marker = TRUE, marker =custom_marker )


# Run scCatch to predict cell types
scCatch_obj <- scCATCH.predict(scCatch_obj)

# Extract predicted cell types
predicted_cell_types <- scCatch_obj@meta.data$cell_type

# Add predicted cell types to the Seurat object metadata
veh.seurat.obj$predicted_cell_type <- predicted_cell_types

# Visualize the annotated results
DimPlot(veh.seurat.obj, reduction = "umap", group.by = "predicted_cell_type")


























sc_data@meta.data


gene_signatures <- apply(gene_signatures, 2, toupper)

print(str(gene_signatures))
print(head(gene_signatures))




##to check whether genes are present or not in both 
# Function to check presence of gene signatures in Seurat object
check_genes_presence <- function(gene_list, seurat_obj) {
  sapply(gene_list, function(genes) all(genes %in% rownames(seurat_obj)))
}

# Apply this function
gene_presence <- check_genes_presence(gene_signatures, veh.seurat.obj)
print(gene_presence)


# Check what's left
print(names(sufficient_features))




# Assuming 'sc_data' is your Seurat object and 'gene_signatures' is a list of your gene sets
sc_data <- AddModuleScore(object = veh.seurat.obj, features = gene_signatures, name = 'signature_score')



view(sc_data)
---------------------------------------------------------------------------
#using SingleR for annotation
ref <- celldex::HumanPrimaryCellAtlasData()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

results <- SingleR(test = as.SingleCellExperiment(veh.seurat.obj), ref = ref, labels = ref$label.main)


veh.seurat.obj$singlr_labels <- results$labels

DimPlot(veh.seurat.obj, reduction = 'umap', group.by = 'singlr_labels', label = TRUE)









































# Subset the metadata to include only cells from cluster 0
cluster_0_metadata <- veh.seurat.obj@meta.data[veh.seurat.obj@meta.data$seurat_clusters == 0, ]

# Print the subsetted metadata
print(cluster_0_metadata)


# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA",only.pos = TRUE)
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_0.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 1, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA",only.pos = TRUE)
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_1.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 2, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA")
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_2.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 3, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA")
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_3.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 4, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA")
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_4.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 5, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA")
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_5.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 6, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA")
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_6.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 7, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA")
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_7.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 8, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA")
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_8.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 9, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA")
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_9.csv", row.names = TRUE)
# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 10, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA")
# Save the 'markers' data frame to a CSV file
write.csv(markers, "markers_10.csv", row.names = TRUE)











# Load necessary library
library(dplyr)

# Read the CSV files
cluster0_genes <- read.csv("markers_9.csv", check.names = FALSE)  # Make sure to include 'check.names = FALSE' if there are special characters in column names

# Assuming 'signature_markers' is a separate CSV file
signature_markers <- read.csv("genesignature_naoto.csv")

# Rename columns for clarity if necessary
colnames(cluster0_genes) <- c("Gene", "myAUC", "avg_diff", "power", "avg_log2FC", "pct.1", "pct.2")

# Merge the data frames based on gene names
merged_data <- merge(cluster0_genes, signature_markers, by.x = "Gene", by.y = "Name")

if (nrow(merged_data) > 0) {
  # Write the merged data to a CSV file
  write.csv(merged_data, "annotated_9.csv", row.names = FALSE)
  print("Merged data has been saved to 'annotated_1.csv'.")
} else {
  print("No matching genes found.")
}

----------------------------------------------------
# Read the CSV files
cluster0_genes <- read.csv("markers_1.csv", check.names = FALSE)  # Ensure 'check.names = FALSE' for special characters in column names

# Assuming 'signature_markers' is a separate CSV file
signature_markers <- read.csv("genesignature_naoto.csv")

# Convert gene names in both data frames to upper case for case-insensitive matching
cluster0_genes$Gene <- toupper(cluster0_genes$Gene)
signature_markers$Name <- toupper(signature_markers$Name)

# Rename columns for clarity if necessary
colnames(cluster0_genes) <- c("Gene", "myAUC", "avg_diff", "power", "avg_log2FC", "pct.1", "pct.2")

# Merge the data frames based on gene names
merged_data <- merge(cluster0_genes, signature_markers, by.x = "Gene", by.y = "Name")

if (nrow(merged_data) > 0) {
  # Write the merged data to a CSV file
  write.csv(merged_data, "annotated_1.csv", row.names = FALSE)
  print("Merged data has been saved to 'annotated_9.csv'.")
} else {
  print("No matching genes found.")
}
-------

  # Load necessary library
  library(dplyr)

# Read the CSV file
data <- read.csv("annotated_1.csv")

# Group by 'list' column and summarize counts
list_counts <- data %>%
  group_by(list) %>%
  summarise(Count = n())

# Print the counts for each unique entry in the 'list' column
print(list_counts)

# Optionally, you can also write this output to a CSV file
write.csv(list_counts, "list_counts.csv", row.names = FALSE)

----------------------------------------
  
# Read the CSV files
cluster0_genes <- read.csv("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/scRNAseq/Veh1/outs/cluster1.csv", check.names = FALSE)  # Ensure 'check.names = FALSE' for special characters in column names

# Assuming 'signature_markers' is a separate CSV file
signature_markers <- read.csv("genesignature_naoto.csv")

# Convert gene names in both data frames to upper case for case-insensitive matching
cluster0_genes$Gene <- toupper(cluster0_genes$FeatureName)
signature_markers$Name <- toupper(signature_markers$Name)

# Rename columns for clarity if necessary
colnames(cluster0_genes) <- c("FeatureName",	"Cluster 1 Average",	"Cluster 1 Log2 Fold Change")

# Merge the data frames based on gene names
merged_data <- merge(cluster0_genes, signature_markers, by.x = "FeatureName", by.y = "Name")

if (nrow(merged_data) > 0) {
  # Write the merged data to a CSV file
  write.csv(merged_data, "annotated_1_loupe.csv", row.names = FALSE)
  print("Merged data has been saved to 'annotated_9.csv'.")
} else {
  print("No matching genes found.")
} 

----------------------------------------------------

# Find markers
markers <- FindMarkers(veh.seurat.obj, ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc", assay = "RNA", only.pos = TRUE)

# Sort the data frame by 'avg_log2FC' in decreasing order and select the top 50 genes
top_markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ][1:50, ]

# Save the 'top_markers' data frame to a CSV file
write.csv(top_markers, "markers_0_TOP50.csv", row.names = TRUE)

# Read the CSV files
cluster0_genes <- read.csv("markers_0_TOP50.csv", check.names = FALSE)  # Ensure 'check.names = FALSE' for special characters in column names

# Assuming 'signature_markers' is a separate CSV file
signature_markers <- read.csv("genesignature_naoto.csv")

# Convert gene names in both data frames to upper case for case-insensitive matching
cluster0_genes$Gene <- toupper(cluster0_genes$Gene)
signature_markers$Name <- toupper(signature_markers$Name)

# Rename columns for clarity if necessary
colnames(cluster0_genes) <- c("Gene", "myAUC", "avg_diff", "power", "avg_log2FC", "pct.1", "pct.2")

# Merge the data frames based on gene names
merged_data <- merge(cluster0_genes, signature_markers, by.x = "Gene", by.y = "Name")

if (nrow(merged_data) > 0) {
  # Write the merged data to a CSV file
  write.csv(merged_data, "annotated_0_top50.csv", row.names = FALSE)
  print("Merged data has been saved to 'annotated_9.csv'.")
} else {
  print("No matching genes found.")
}

data <- read.csv("annotated_0_top50.csv")

# Group by 'list' column and summarize counts
list_counts <- data %>%
  group_by(list) %>%
  summarise(Count = n())

# Print the counts for each unique entry in the 'list' column
print(list_counts)

# Optionally, you can also write this output to a CSV file
write.csv(list_counts, "list_counts.csv", row.names = FALSE)
