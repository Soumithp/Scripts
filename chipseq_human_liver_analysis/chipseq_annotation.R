# ChIP-seq peak annotation + functional enrichment in R (ChIPseeker / clusterProfiler).
# Companion to chipseq_analysis.ipynb — same H3K27ac human-liver peaks, R side of the workflow.
# Input: peaks.bed (ENCODE narrowPeak, GRCh38) downloaded by the notebook.

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile("peaks.bed")

# annotate peaks relative to gene features (promoter, UTR, intron, distal, ...)
peakAnno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000),
                         annoDb = "org.Hs.eg.db")

dir.create("results", showWarnings = FALSE)

# feature distribution + distance-to-TSS
png("results/R_peak_annotation_pie.png", width = 900, height = 700, res = 130)
plotAnnoPie(peakAnno); dev.off()

png("results/R_dist_to_TSS.png", width = 1000, height = 500, res = 130)
plotDistToTSS(peakAnno); dev.off()

# GO enrichment on genes near peaks
genes <- unique(as.data.frame(peakAnno)$geneId)
ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "BP",
                pvalueCutoff = 0.05, readable = TRUE)

png("results/R_GO_enrichment.png", width = 1100, height = 700, res = 130)
dotplot(ego, showCategory = 15) + ggtitle("GO (BP) — genes near H3K27ac peaks")
dev.off()

write.csv(as.data.frame(peakAnno), "results/R_peak_annotation.csv", row.names = FALSE)
write.csv(as.data.frame(ego),      "results/R_GO_enrichment.csv",   row.names = FALSE)
cat("done — see results/\n")
