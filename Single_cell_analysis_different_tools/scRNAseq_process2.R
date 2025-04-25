# Required libraries
library(Seurat)
library(dplyr)
library(openxlsx)
library(cowplot)

setwd("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki")
# Step 1: Read and preprocess the Seurat object
svr_data <- Read10X_h5(filename = "C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/hg38/hg38_filtered_feature_bc_matrix.h5")
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


# Step 1: Read Gene Signatures
gene_signatures <- read.table("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiver_MP_signature.txt", header = TRUE, sep = "\t", check.names = FALSE)

# Convert gene signatures to a list of gene sets
gene_list <- lapply(gene_signatures, as.character)
features_list <- lapply(colnames(gene_signatures), function(x) gene_signatures[[x]])
names(features_list) <- colnames(gene_signatures)

































# Step 1: Read Gene Signatures
gene_signatures <- read.table("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/scLiverMP_Tcell_signature.txt", header = TRUE, sep = "\t", check.names = FALSE)

# Convert gene signatures to a list of gene sets
gene_list <- lapply(gene_signatures, as.character)
features_list <- lapply(colnames(gene_signatures), function(x) gene_signatures[[x]])
names(features_list) <- colnames(gene_signatures)

# Prepend "A_" to each feature name
#features_list <- setNames(features_list, paste0("A_", names(features_list)))
# Step 2: Calculate AddModuleScore with custom names
# Step 2: Calculate AddModuleScore
svr.seurat.obj <- AddModuleScore(svr.seurat.obj, features = features_list, name = "ModuleScore")

# Step 3: Extract Module Scores
#module_scores <- svr.seurat.obj@meta.data %>%dplyr::select(starts_with("A_"))
# Step 4: Save the Module Scores to a text file
#write.table(module_scores, file = "module_scores_tcell.txt", sep = "\t", quote = FALSE, col.names = NA)


# Create a list to store the UMAP plots
umap_plots <- list()

#Visualize module scores for each T-cell gene signature and store plots
for (i in 1:length(features_list)) {
  feature_name <- paste0("ModuleScore", i)
  p <- FeaturePlot(svr.seurat.obj, features = feature_name) +
    ggplot2::ggtitle(paste("UMAP plot for", names(features_list)[i])) +
    NoLegend()
  umap_plots[[i]] <- p
}
# Split the UMAP plots into two groups
umap_plots_1 <- umap_plots[1:15]
umap_plots_2 <- umap_plots[26:length(umap_plots)]

# Combine the first group of UMAP plots into one frame
combined_umap_plot_1 <- plot_grid(plotlist = umap_plots_1, ncol = 5, align = 'hv')

# Combine the second group of UMAP plots into one frame
combined_umap_plot_2 <- plot_grid(plotlist = umap_plots_2, ncol = 5, align = 'hv')

# Print the combined UMAP plots
print(combined_umap_plot_1)
print(combined_umap_plot_2)

# Save the combined UMAP plots
ggsave(filename = "module_scores_tcellsig_umap_plot_1.png", plot = combined_umap_plot_1, width = 16, height = 12)
ggsave(filename = "module_scores_tcellsigumap_plot_2.png", plot = combined_umap_plot_2, width = 16, height = 12)



# Create a list to store the violin plots
vln_plots <- list()

# Create violin plots for module scores by cluster and store plots
for (i in 1:length(features_list)) {
  feature_name <- paste0("ModuleScore", i)
  p <- VlnPlot(svr.seurat.obj, features = feature_name) +
    ggplot2::ggtitle(paste("for", names(features_list)[i])) +
    theme(plot.title = element_text(size = 10))  # Adjust title size
  vln_plots[[i]] <- p
}

# Split the violin plots into two groups
vln_plots_1 <- vln_plots[1:15]
vln_plots_2 <- vln_plots[26:length(vln_plots)]

# Combine the first group of violin plots into one frame
combined_vln_plot_1 <- plot_grid(plotlist = vln_plots_1, ncol = 5, align = 'hv')

# Combine the second group of violin plots into one frame
combined_vln_plot_2 <- plot_grid(plotlist = vln_plots_2, ncol = 5, align = 'hv')

# Print the combined violin plots
print(combined_vln_plot_1)
print(combined_vln_plot_2)

# Save the combined violin plots
ggsave(filename = "module_scores_tcell_vln_plot_1.png", plot = combined_vln_plot_1, width = 16, height = 12)
ggsave(filename = "module_scores_vln_plot_2.png", plot = combined_vln_plot_2, width = 16, height = 12)



# Summarize module scores by cluster
module_scores_by_cluster <- svr.seurat.obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(starts_with("ModuleScore"), mean))

write.table(module_scores_by_cluster, file = "module_scores_bycluster(mean).txt", sep = "\t", quote = FALSE, col.names = NA)

# Print the summary
print(module_scores_by_cluster)


# Step 4: Find the highest module score for each cluster
highest_scores <- module_scores_by_cluster %>%
  group_by(seurat_clusters) %>%
  summarize(across(starts_with("ModuleScore"), max, .names = "max_{col}")) %>%
  pivot_longer(cols = starts_with("max_ModuleScore"), names_to = "ModuleScore", values_to = "Score") %>%
  group_by(seurat_clusters) %>%
  slice_max(Score) %>%
  ungroup()

# Step 6: Annotate clusters with the highest scoring module
cluster_annotations <- highest_scores %>%
  dplyr::select(seurat_clusters, ModuleScore)

# View the cluster annotations
print(cluster_annotations)






# Add the annotations to the Seurat object's metadata
svr.seurat.obj@meta.data <- svr.seurat.obj@meta.data %>%
  left_join(cluster_annotations, by = "seurat_clusters")
# Create a named vector of cluster annotations based on the provided data
cluster_labels <- c(
  "0" = "HepatoImmunoactive",
  "1" = "HepatoImmunoactive",
  "2" = "HepatoImmunoactive",
  "3" = "HepatoImmunoactive",
  "4" = "Cholangiocyte",
  "5" = "Endoth_CV_LSEC"
)

# Create a new column in the metadata for cluster annotations
svr.seurat.obj@meta.data$ClusterAnnotation <- as.character(svr.seurat.obj@meta.data$seurat_clusters)

# Assign the annotations based on the cluster IDs
svr.seurat.obj@meta.data$ClusterAnnotation <- recode(svr.seurat.obj@meta.data$ClusterAnnotation, !!!cluster_labels)

# Convert ClusterAnnotation to a factor
svr.seurat.obj@meta.data$ClusterAnnotation <- as.factor(svr.seurat.obj@meta.data$ClusterAnnotation)

# Print the metadata to verify the annotations
print(head(svr.seurat.obj@meta.data, 20))

# Check for NA values in the ClusterAnnotation column
if (any(is.na(svr.seurat.obj@meta.data$ClusterAnnotation))) {
  stop("Some cells have not been assigned a cluster annotation.")
}

# Adjust scaling parameters (adjust limits as needed)
p <- DimPlot(svr.seurat.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = "ClusterAnnotation",pt.size = 0.5) + ggtitle("UMAP by Cluster Annotation")

p <- p + scale_x_continuous(limits = c(-7, 7))  # Adjust x-axis limits
p <- p + scale_y_continuous(limits = c(-7, 7))  # Adjust y-axis limits

# Print the UMAP plot
print(p)












# Step 3: Identify the highest scoring gene signature for each cluster
module_scores_long <- module_scores_by_cluster %>%
  pivot_longer(cols = starts_with("ModuleScore"), names_to = "signature", values_to = "score")

# Extract the highest scoring signature for each cluster
highest_scoring_signatures <- module_scores_long %>%
  group_by(seurat_clusters) %>%
  slice_max(order_by = score, n = 1) %>%
  ungroup() %>%
  mutate(signature_name = sub("ModuleScore", "", signature))

# Create a mapping from cluster to highest scoring gene signature name
cluster_to_signature <- highest_scoring_signatures %>%
  select(seurat_clusters, signature_name)

# Step 4: Annotate clusters based on highest module score
svr.seurat.obj$cluster_annotation <- svr.seurat.obj$seurat_clusters
svr.seurat.obj$cluster_annotation <- as.character(svr.seurat.obj$cluster_annotation)

for (i in 1:nrow(cluster_to_signature)) {
  cluster <- cluster_to_signature$seurat_clusters[i]
  signature_name <- cluster_to_signature$signature_name[i]
  svr.seurat.obj$cluster_annotation[svr.seurat.obj$seurat_clusters == cluster] <- signature_name
}

# Step 5: Visualize the annotated clusters
DimPlot(svr.seurat.obj, group.by = "cluster_annotation", label = TRUE, repel = TRUE) +
  ggplot2::ggtitle("Annotated Clusters by Gene Signature") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Save the annotated UMAP plot
ggsave(filename = "annotated_umap_plot.png", width = 10, height = 8)

















###########just to visulaize ddr1
# Gene of interest
gene_of_interest <- "DDR"

# Step 1: Feature Plot for UMAP
feature_plot <- FeaturePlot(svr.seurat.obj, features = gene_of_interest, reduction = "umap") +
  ggplot2::ggtitle(paste("UMAP plot for", gene_of_interest))
print(feature_plot)

# Step 2: Heatmap
# Retrieve the expression data for the gene of interest
expression_data <- FetchData(svr.seurat.obj, vars = gene_of_interest)

# Add cluster information to the expression data
expression_data$Cluster <- Idents(svr.seurat.obj)

# Create a heatmap for the gene of interest
heatmap_plot <- DoHeatmap(svr.seurat.obj, features = gene_of_interest) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::ggtitle(paste("Heatmap for", gene_of_interest))
print(heatmap_plot)

























#Find marker genes for each cluster
markers_seur <- FindAllMarkers(svr.seurat.obj, only.pos = TRUE)

# View the markers
head(markers_seur)

# Save the markers to an Excel file
write.xlsx(markers_seur, "cluster_markers.xlsx", row.names = FALSE)


# Step 4: Retrieve the top 5 marker genes per cluster
top5 <- markers_seur %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_seur), value = TRUE)),
                   n = 10)

# Step 5: Create the dot plot
DotPlot(svr.seurat.obj, features = unique(top5$gene)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,
                                                     size = 8, hjust = 1)) +
  NoLegend()

# Step 6: Create the heatmap
DoHeatmap(svr.seurat.obj, features = unique(top5$gene)) +
  NoLegend() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))


##################common genes from gene signature
# Step 4: Read the gene signatures and convert to a list of gene sets
gene_signatures <- read.table("C:/Users/S226953/OneDrive - University of Texas Southwestern/Desktop/ev59_scRNAseq_hiroaki/DDR1hepg2_FDR025_sig_20231018_FINAL.txt", header = TRUE, sep = "\t", check.names = FALSE)
gene_list <- lapply(gene_signatures, as.character)
features_list <- lapply(colnames(gene_signatures), function(x) gene_signatures[[x]])
names(features_list) <- colnames(gene_signatures)

# Flatten the gene signatures to get a unique list of all signature genes
all_signature_genes <- unique(unlist(features_list))

# Step 5: Find the intersection of marker genes and gene signature genes
common_genes <- markers_seur %>% filter(gene %in% all_signature_genes)

# Save the common genes to an Excel file
write.xlsx(common_genes, "common_marker_genes.xlsx", row.names = FALSE)

# Print the common genes
print(common_genes)

# Optionally, create a dot plot and heatmap for the common genes
if (nrow(common_genes) > 0) {
  common_gene_list <- unique(common_genes$gene)
  
  # Dot plot
  DotPlot(svr.seurat.obj, features = common_gene_list) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
    NoLegend()
  
  # Heatmap
  DoHeatmap(svr.seurat.obj, features = common_gene_list) +
    NoLegend() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}


# Find the intersection of marker genes and gene signature genes
common_genes <- markers_seur %>% filter(gene %in% all_signature_genes)

# Retrieve the top 20 common marker genes per cluster based on the average log fold change
top20_common_genes <- common_genes %>% group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 20)

# Save the top 20 common genes to an Excel file
write.xlsx(top20_common_genes, "top20_common_marker_genes.xlsx", row.names = FALSE)

# Print the top 20 common genes
print(top20_common_genes)

# Optionally, create a dot plot and heatmap for the top 20 common genes
if (nrow(top20_common_genes) > 0) {
  top20_common_gene_list <- unique(top20_common_genes$gene)
  
  # Dot plot
  DotPlot(svr.seurat.obj, features = top20_common_gene_list) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
    NoLegend()
  
  # Heatmap
  DoHeatmap(svr.seurat.obj, features = top20_common_gene_list) +
    NoLegend() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}



###########DEGS
# Step 2: Find DEGs for each cluster
deg_list <- list()
clusters <- unique(Idents(svr.seurat.obj))
for (cluster in clusters) {
  degs <- FindMarkers(svr.seurat.obj, ident.1 = cluster, only.pos = TRUE)
  degs <- degs %>% mutate(gene = rownames(degs)) # Add gene column
  degs$cluster <- cluster
  deg_list[[as.character(cluster)]] <- degs
}

# Combine all DEGs into a single data frame
deg_df <- do.call(rbind, lapply(names(deg_list), function(x) {
  deg_list[[x]] %>% mutate(cluster = x)
}))

# Save the DEGs to an Excel file
write.xlsx(deg_df, "deg_per_cluster.xlsx", row.names = FALSE)

# Print the first few rows of the DEG data frame
head(deg_df)

# Step 3: Select the top 20 DEGs for each cluster
top_deg_per_cluster <- deg_df %>% group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 20)

# Save the top DEGs to an Excel file
write.xlsx(top_deg_per_cluster, "top_deg_per_cluster.xlsx", row.names = FALSE)

# Print the top DEGs
print(top_deg_per_cluster)

# Step 4: Create plots for the top DEGs
if (nrow(top_deg_per_cluster) > 0) {
  top_deg_gene_list <- unique(top_deg_per_cluster$gene)
  
  # Dot plot
  dot_plot <- DotPlot(svr.seurat.obj, features = top_deg_gene_list) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, size = 7, hjust = 1)) +
    NoLegend()
  print(dot_plot)
  
  # Heatmap
  heatmap_plot <- DoHeatmap(svr.seurat.obj, features = top_deg_gene_list) +
    NoLegend() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  print(heatmap_plot)
}