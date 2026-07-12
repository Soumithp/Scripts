# Bioinformatics scripts

A collection of the analysis code I've written and used across my computational biology work, mostly in
liver cancer and hepatocellular carcinoma (HCC) genomics. It spans bulk RNA-seq, single-cell and
single-nucleus RNA-seq, spatial transcriptomics, ChIP-seq, viral integration, and somatic mutation
analysis — plus survival/prognostic modeling and a few self-contained, reproducible analysis projects.

The repo is split into two kinds of thing: **analysis projects** that run end-to-end and ship with
their own results and write-ups, and **pipelines & scripts** — the working code behind larger studies.

## Analysis projects (notebooks + results)

Self-contained, documented, and reproducible — each runs top-to-bottom on Google Colab and has figures
and a walkthrough in its own README.

| Project | What it shows | Tools |
|---|---|---|
| [`scrna_celltype_ML_classification`](scrna_celltype_ML_classification) | scRNA-seq clustering + ML cell-type classification on ~11k PBMCs (~99% accuracy) | Scanpy, scikit-learn |
| [`chipseq_human_liver_analysis`](chipseq_human_liver_analysis) | H3K27ac active-regulatory landscape in primary human liver (ENCODE) — TSS heatmaps, peak annotation, GO | deepTools, pybedtools, ChIPseeker |
| [`visium_spatial_analysis`](visium_spatial_analysis) | 10x Visium spatial domains, marker genes, spatially variable genes, neighborhood structure | Scanpy, Squidpy |

## Pipelines & scripts

Working code from research projects. These generally have hardcoded cluster paths — update paths/inputs
before reuse. Each folder has a README with the run order.

| Folder | What it does | Main tools |
|---|---|---|
| [`RNAseq_expression`](RNAseq_expression) | Bulk RNA-seq end-to-end: QC → alignment → counting → filtering → normalization → differential expression | STAR, HISAT2, featureCounts, DESeq2, edgeR |
| [`Single_cell_analysis_different_tools`](Single_cell_analysis_different_tools) | scRNA-seq / snRNA-seq processing and cell-type annotation compared across tools | Seurat, SingleR, scType, CellAssign, Azimuth |
| [`spatial_transcriptomics`](spatial_transcriptomics) | 10x Visium preprocessing (FASTQ generation + FFPE image-to-count) | Space Ranger |
| [`HBV_fusion`](HBV_fusion) | Detect HBV viral integration / host–virus fusion events in RNA-seq | Kraken, ViFi / FastViFi |
| [`HBV_signature`](HBV_signature) | Prognostic gene signature for HBV-related HCC, validated with survival analysis | NTP, LOOCV, Cox, Kaplan-Meier |
| [`SNP_mutational_patterns`](SNP_mutational_patterns) | Somatic/germline variant calling → COSMIC mutational signatures | GATK, Picard, ANNOVAR, MutationalPatterns |

## Languages & tools

R (Seurat, DESeq2, edgeR, WGCNA, MutationalPatterns, ChIPseeker, survival), Python (pandas, numpy,
pysam, scanpy, squidpy, scikit-learn, deepTools), and bash for the HPC/pipeline steps (SLURM modules,
STAR/HISAT2, GATK, Space Ranger, FastViFi). ANNOVAR (Perl) is used for variant annotation.

## A note on third-party code

Some scripts wrap or include tools written by others; I've kept their original author credits rather
than relabeling them:

- **ANNOVAR** (`prepare_annovar_user.pl`) — Kai Wang
- **Nearest Template Prediction** (`NTP.R`) — Yujin Hoshida, Broad Institute
- **`strat.split.stef.R`** — original stratified-split utility ("Stef", 2007)

Everything else is my own analysis code from these projects.
