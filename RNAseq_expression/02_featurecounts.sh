#!/bin/bash
# featureCounts: BAMs -> gene count matrix (feeds RawCountSummary.R downstream)
module load subread/2.0.1

threads=8
gtf=/path/to/annotation.gtf
bamDir=/path/to/aligned
outDir=/path/to/counts
mkdir -p "$outDir"

# paired-end (-p --countReadPairs), count primary alignments only, gene level
featureCounts -T "$threads" -p --countReadPairs --primary -t exon -g gene_id \
    -a "$gtf" \
    -o "$outDir/featurecounts.primary.txt" \
    "$bamDir"/*Aligned.sortedByCoord.out.bam

echo "counts -> $outDir/featurecounts.primary.txt"
# RawCountSummary.R / RawCountSummary_batchcode.R takes this featurecounts.primary.txt as input.
