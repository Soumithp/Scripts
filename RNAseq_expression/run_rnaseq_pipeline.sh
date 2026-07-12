#!/bin/bash
# End-to-end RNA-seq expression pipeline.
# QC -> alignment -> featureCounts -> count summary -> filtering -> normalization -> DEG.
# Set the variables below, then run: bash run_rnaseq_pipeline.sh
set -euo pipefail

# ---- project settings ----
projectName="my_project"
fastqDir=/path/to/fastq
gtf=/path/to/annotation.gtf
annotTable=/path/to/gene_annotation.txt      # gene_id -> symbol/length table used by the R steps
workDir=/path/to/work
countsDir="$workDir/counts"
scriptDir=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$workDir"

# ---- 1. read QC ----
bash "$scriptDir/fastqc.sh"
bash "$scriptDir/mutliqc.sh"

# ---- 2. alignment (STAR primary; swap to 01b_align_hisat2.sh for HISAT2) ----
bash "$scriptDir/01_align_star.sh"

# ---- 3. featureCounts ----
bash "$scriptDir/02_featurecounts.sh"

# ---- 4. raw count summary (combines per-sample featureCounts into one matrix) ----
Rscript --vanilla "$scriptDir/RawCountSummary_batchcode.R" \
    "$scriptDir/RawCountSummary.R" "$projectName" "$annotTable" \
    "$countsDir/featurecounts.primary.txt"

# ---- 5. filtering: drop mostly-zero genes, then low-variance genes ----
Rscript --vanilla "$scriptDir/ZeroProp_batch.R"  "$scriptDir/ZeroProp.R"  "$projectName"
Rscript --vanilla "$scriptDir/VarFilter_batch.R" "$scriptDir/VarFilter.R" "$projectName"

# ---- 6. normalization ----
Rscript --vanilla "$scriptDir/normalizeData_batch.R" "$scriptDir/normalizeData.R" "$projectName"

# ---- 7. differential expression (DESeq2 / edgeR) ----
Rscript --vanilla "$scriptDir/DEG_analysis.R"

echo "RNA-seq pipeline finished for $projectName"
# Note: the R step arguments follow each *_batch.R script's own usage; adjust names/paths to your run.
