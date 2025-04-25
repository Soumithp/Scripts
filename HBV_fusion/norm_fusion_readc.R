library(DESeq2)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3 ) {
  stop("Usage: Rscript --vanilla norm_fusion_readc.R <CountsAll_gct> <CountsFusion_csv> <output>", call.=FALSE)
} else if (length(args) == 3) {
  inputFile <- args[1]
  inputHBV <- args[2]
  outFile <- args[3]
}

# whole transcriptome matrix
# inputFile <- "/project/shared/DSSR/s229294/hoshida_lab_git/hbv_fusion_detection/workflow/output/merged_Rawcountsmatrix.gct"
countMat <- read.delim(inputFile, header = TRUE, skip = 2, sep = "\t")
rownames(countMat) <- countMat[,1]
countMat <- countMat[,-c(1,2)]
# remove samples with 0 counts for all genes
countMat <- countMat[,colnames(countMat)[colSums(countMat) > 0]]

# hbv fusion count matrix
# inputHBV <- "/project/shared/DSSR/s229294/hoshida_lab_git/hbv_fusion_detection/workflow/output/hbvfusion_summary_readc_matrix.csv"
hbvMat <- read.delim(inputHBV, header = TRUE, sep = ",")
rownames(hbvMat) <- hbvMat[,1]
hbvMat <- hbvMat[,-1]

# generate size factor
sizeFactor <- estimateSizeFactorsForMatrix(countMat)
sizeFactor <- sizeFactor[colnames(hbvMat)]

norm_hbvMat <- sweep(hbvMat, 2, sizeFactor, FUN = "/")
norm_hbvMat <- round(norm_hbvMat, digits = 2)
# outFile <- "/project/shared/DSSR/s229294/hoshida_lab_git/hbv_fusion_detection/workflow/output/hbvfusion_summary_readc_matrix_norm.csv"
write.csv(norm_hbvMat, file = outFile,
          quote = F)
