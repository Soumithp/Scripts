# Bioinformatics scripts

This is a collection of the analysis code I've written and used across my computational biology work,
mostly in liver cancer and hepatocellular carcinoma (HCC) genomics. It spans several data types —
bulk RNA-seq, single-cell and single-nucleus RNA-seq, spatial transcriptomics, viral integration, and
somatic mutation analysis — plus survival/prognostic modeling and a small end-to-end machine-learning
project.

These are real working scripts, not a packaged tool. Many of them have file paths from the cluster or
machine I ran them on, so if you want to reuse something, expect to update the paths and inputs first.
Each folder has its own README with the details.

## What's here

| Folder | What it does | Main tools |
|---|---|---|
| [`RNAseq_expression`](RNAseq_expression) | Bulk RNA-seq: QC → counting → filtering → normalization → differential expression | FastQC, MultiQC, DESeq2, edgeR |
| [`Single_cell_analysis_different_tools`](Single_cell_analysis_different_tools) | scRNA-seq / snRNA-seq processing and cell-type annotation compared across tools | Seurat, SingleR, scType, CellAssign, Azimuth |
| [`spatial_transcriptomics`](spatial_transcriptomics) | 10x Visium preprocessing (FASTQ + FFPE image-to-count) | Space Ranger |
| [`HBV_fusion`](HBV_fusion) | Detect HBV viral integration / host–virus fusion events in RNA-seq | Kraken, ViFi / FastViFi |
| [`HBV_signature`](HBV_signature) | Prognostic gene signature for HBV-related HCC, validated with survival analysis | NTP, LOOCV, Cox, Kaplan-Meier |
| [`SNP_mutational_patterns`](SNP_mutational_patterns) | Somatic variants → COSMIC mutational signatures | MutationalPatterns, ANNOVAR |
| [`scrna_celltype_ML_classification`](scrna_celltype_ML_classification) | End-to-end scRNA-seq clustering + ML cell-type classification (documented, reproducible) | Scanpy, scikit-learn |

If you're just browsing, [`scrna_celltype_ML_classification`](scrna_celltype_ML_classification) is the
most self-contained — it runs top to bottom on Colab and has its own results and write-up.

## Languages & tools

R (Seurat, DESeq2, edgeR, WGCNA, MutationalPatterns, survival), Python (pandas, numpy, pysam, scanpy,
scikit-learn), and bash for the HPC/pipeline steps (SLURM modules, Space Ranger, FastViFi). ANNOVAR
(Perl) is used for variant annotation.

## A note on third-party code

Some scripts wrap or include tools written by others, and I've kept their original author credits
rather than relabeling them:

- **ANNOVAR** (`prepare_annovar_user.pl`) — Kai Wang
- **Nearest Template Prediction** (`NTP.R`) — Yujin Hoshida, Broad Institute
- **`strat.split.stef.R`** — original stratified-split utility ("Stef", 2007)

Everything else is my own analysis code from these projects.
