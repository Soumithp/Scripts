# Bulk RNA-seq expression analysis

An end-to-end bulk RNA-seq workflow — from raw FASTQ through alignment, counting, filtering,
normalization, and differential expression (DESeq2 / edgeR). `run_rnaseq_pipeline.sh` chains the
whole thing; you can also run any step on its own. Most R steps have a base script plus a `_batch`
version that runs the same step across several datasets.

### End-to-end order

1. `fastqc.sh`, `mutliqc.sh` — read QC (FastQC → MultiQC).
2. `01_align_star.sh` — alignment with STAR (primary). `01b_align_hisat2.sh` is the HISAT2 alternative.
3. `02_featurecounts.sh` — featureCounts: BAMs → `featurecounts.primary.txt`.
4. `RawCountSummary.R` (`_batchcode`) — combine per-sample counts into one matrix + summary stats.
5. `ZeroProp.R` (`_batch`) — drop genes with too many zeros.
6. `VarFilter.R` (`_batch`) — variance-based gene filtering.
7. `normalizeData.R` (`_batch`) — normalization (with `logFunction.R` helper).
8. `DEG_analysis.R` + `DEG_analysis_function.R` — differential expression.

`IntersampleCorrelation.R` (`_batch`) and `Gene_couting.R` are extra QC/utility steps you can run
alongside (sample correlation check, non-zero gene counts).

### Run it

```bash
# edit the variables at the top first (paths, GTF, project name)
bash run_rnaseq_pipeline.sh
```

### Notes

- The alignment/count scripts use `module load` and reference paths — set them for your cluster.
- The `_batch` wrappers `source()` the base scripts, so keep each pair together.
