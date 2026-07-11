# scRNA-seq PBMC Cell-Type Classification

## Abstract

This project builds a complete single-cell RNA-seq (scRNA-seq) analysis pipeline on ~11,000 human
peripheral blood mononuclear cells (PBMCs). Cells are quality-controlled, normalized, clustered
(Leiden) and visualized (UMAP); each cluster is annotated to an immune lineage using curated marker
panels of 15+ genes, and every panel is confirmed statistically enriched in its cluster (Wilcoxon
rank-sum, Benjamini–Hochberg adjusted p < 0.05). A supervised machine-learning step then trains a
Random Forest and an MLP to predict cell type directly from gene expression, reaching **98.7% test
accuracy** (macro-F1 0.98) with **5-fold cross-validation of 0.984 ± 0.004**. Top model features
recover the canonical marker genes, confirming the classifiers learned real biology rather than
technical artifacts.

## Data

| | |
|---|---|
| **Dataset** | 10k PBMCs from a Healthy Donor, v3 chemistry (10x Genomics, public) |
| **Assay origin** | Droplet scRNA-seq (Zheng et al., *Nature Communications* 2017) |
| **Cells analyzed** | 10,685 (after QC) |
| **Genes** | 20,292 (after gene filtering) |
| **Access** | Downloaded in the notebook, no login |

> Zheng GXY, Terry JM, Belgrader P, et al. *Massively parallel digital transcriptional profiling of
> single cells.* **Nature Communications** 8, 14049 (2017). https://doi.org/10.1038/ncomms14049

## Pipeline & parameters

| Step | Method | Key parameters |
|---|---|---|
| Quality control | filter cells/genes, drop doublets & high-mito | `min_genes=200`, `min_cells=3`, `n_genes_by_counts < 6000`, `pct_counts_mt < 15%` |
| Normalization | counts-per-10k + log1p | `target_sum=1e4`, then `log1p` |
| Feature selection | highly variable genes | `n_top_genes=2000`, `flavor='seurat'` |
| Scaling | z-score, clipped | `max_value=10` |
| Dimensionality reduction | PCA | `n_comps=50`, `svd_solver='arpack'`, 30 PCs used downstream |
| Graph & clustering | kNN graph + Leiden | `n_neighbors=15`, `n_pcs=30`, `resolution=0.5` → 19 clusters |
| Visualization | UMAP | default, `random_state=0` |
| Marker detection | Wilcoxon rank-sum | `method='wilcoxon'`, BH-adjusted p |
| Annotation | marker-panel scoring | `sc.tl.score_genes` per panel, argmax label |
| Classification | Random Forest / MLP | RF `n_estimators=300`; MLP `hidden_layer_sizes=(128,64)`, `early_stopping=True` |
| Validation | held-out test + CV | stratified 75/25 split, 5-fold stratified CV |

## Selection criteria & statistics

- **Cell QC:** cells kept if they express ≥200 genes and <6,000 genes (upper bound removes probable
  doublets), and mitochondrial content <15% (high mito fraction marks stressed/lysing cells). Genes
  kept if detected in ≥3 cells.
- **Feature selection:** the top 2,000 highly variable genes drive clustering and are the ML feature
  set — this focuses the model on informative genes and reduces noise/dimensionality.
- **Marker significance:** for each annotated cell type, its reference markers are tested with a
  Wilcoxon rank-sum test (that type vs. all others) and corrected for multiple testing with
  Benjamini–Hochberg. Only markers with **adjusted p < 0.05** are reported. In this run every lineage
  retained **15–20 significant markers**, with adjusted p-values effectively 0 (max ≈ 7.5e-6). Full
  table: [`results/marker_significance.csv`](results/marker_significance.csv).

## Cell types & representative markers

Eight immune populations were resolved. Top significant markers (by adjusted p-value):

| Cell type | Significant markers | Top markers |
|---|---|---|
| B | 19 | CD79A, MS4A1, CD79B, CD37, CD74 |
| CD14+ Monocytes | 20 | S100A9, S100A8, LYZ, MNDA, VCAN |
| CD4 T | 17 | LDHB, TRAC, CD3D, IL7R, TCF7 |
| CD8 T | 15 | NKG7, GZMA, CCL5, KLRG1, IL32 |
| Dendritic cells | 20 | HLA-DPB1, HLA-DRA, HLA-DRB1, PLD4, CST3 |
| FCGR3A+ Monocytes | 19 | LST1, FCGR3A, AIF1, COTL1, FCER1G |
| NK | 19 | NKG7, GNLY, PRF1, KLRD1, CTSW |
| Platelet | 19 | PF4, PPBP, NRGN, CAVIN2, TUBB1 |

Reference panels curated from **PanglaoDB**, **CellMarker 2.0**, and the **Azimuth PBMC** reference.

## Results

| Model | Accuracy | Macro-F1 |
|---|---|---|
| Random Forest | 0.987 | 0.977 |
| MLP | 0.986 | 0.980 |
| Random Forest (5-fold CV) | 0.984 ± 0.004 | — |

Per-class metrics: [`results/classification_report.txt`](results/classification_report.txt) — every
cell type scores F1 ≥ 0.95. Figures and a fuller walk-through are in [`RESULTS.md`](RESULTS.md).

## Interview notes — key decisions & rationale

- **Why cluster first, then classify?** Cell types aren't labeled in raw data. We derive labels by
  unsupervised clustering + marker biology, then train a classifier that could label new cells fast
  without re-clustering.
- **Why a marker *reference*?** The algorithm can't know which genes mark which cell type — we supply
  curated panels (PanglaoDB / CellMarker / Azimuth). The Wilcoxon + BH step is the evidence that
  those panels are actually enriched here.
- **Why macro-F1, not just accuracy?** Classes are imbalanced (CD4 T ~ thousands, platelets ~ tens).
  Macro-F1 weights every class equally, so it catches a model that ignores rare types.
- **How do we know it's not overfitting?** Held-out test set + 5-fold CV (tight ±0.004 spread) +
  feature importances that match known markers (CD3D, MS4A1, NKG7, LYZ, FCGR3A…).
- **Why 15% mito / 6,000-gene cutoffs?** Standard 10x PBMC QC — high mito = dying cells, very high
  gene counts = likely doublets.

## How to run

1. Open `scrna_analysis.ipynb` in Google Colab.
2. Run the install cell, then **Runtime → Restart session**.
3. **Run all**. Figures and tables are written to `results/`.

## Files

- `scrna_analysis.ipynb` — full analysis notebook
- `RESULTS.md` — figures + metrics
- `requirements.txt` — dependencies
- `results/` — figures, `marker_significance.csv`, `classification_report.txt`, `metrics.txt`
