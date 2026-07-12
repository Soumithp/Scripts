# Single-cell / single-nucleus RNA-seq — annotation with different tools

A collection of scRNA-seq (and snRNA-seq) scripts built around Seurat, with cell-type annotation
tried several different ways so I could compare them: SingleR, scType, CellAssign, and Azimuth. Most
start from a 10x matrix/H5, do the usual QC and clustering, then branch off into whichever annotation
approach.

### Preprocessing

- `scRNAseq_from_matrix.R` — read a 10x matrix/H5 into Seurat, QC (% mito), normalize, cluster.
- `scRNAseq_process2.R` — an alternate processing pass.
- `scRNAseq_from_sunny.R` — processing for a collaborator's dataset.
- `snrnaseq.R` — single-nucleus variant of the workflow.

### Annotation approaches

- `scRNAseq_singleR_annot.R`, `scRNAseq_Epacad_singleR_annot.R` — SingleR reference-based annotation.
- `scRNAseq_sctype_annot.R` — scType marker-based annotation.
- `scRNAseq_annot_cellassign.R`, `cellassign.R`, `inference-tensorflow.R`, `utils.R` — CellAssign
  (probabilistic, marker-based; uses TensorFlow).
- `azimuth_analysis.R`, `azimuth_analysis_v2.R` — Azimuth reference mapping (two versions with
  different QC thresholds).

### Signatures, scoring, and utilities

- `changing_mouse_genes_2humangenes.R` — mouse → human gene-symbol conversion.
- `making_list_of_gene_columns.R`, `merging_scores_all&hepato_barcode_based.R` — build/merge
  per-barcode gene-signature scores.
- `scRNAseq_genes_comparing_w_genesig.R` — compare cluster genes against a gene signature.
- `scRNAseq_vis_heatmap.R` — heatmap visualization.
- `simulate.R` — simulated data for testing.

### Notes

- Paths are from my runs; update `setwd()` and file paths first. scType and some references are
  pulled from their public sources at runtime.
