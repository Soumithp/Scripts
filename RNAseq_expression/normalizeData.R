library("DESeq2")
library("edgeR")
normalizeData <- function(inputFile, normalizationMethod, countThres = 0) {
  
  Rawcountsmatrix <- read.delim(inputFile, header = TRUE, skip = 2)
  # Store the first two columns
  gene_info <- Rawcountsmatrix[, 1:2]
  fileName <- tools::file_path_sans_ext(inputFile)
  filename <- file.path(paste0(fileName, "_", normalizationMethod, ".gct"))
  row.names(Rawcountsmatrix)<- Rawcountsmatrix[,1]
  Rawcountsmatrix<- as.matrix(Rawcountsmatrix[, -c(1, 2)])
  
  
  if (normalizationMethod == "RLE") {
    countMatrix <- Rawcountsmatrix[rowSums(Rawcountsmatrix) >= countThres, ]  # No need to filter columns
    colData <- data.frame(rep("null", ncol(countMatrix)))
    dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = colData, design = as.formula('~1'))
    dds <- DESeq(dds)
    normalizedCount <- counts(dds, normalized = TRUE)
  } else if (normalizationMethod == "TMM") {
    countMatrix <- Rawcountsmatrix  # No need to filter columns
    y <- DGEList(countMatrix)
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    normalizedCount <- y$pseudo.counts
  } else {
    stop("Invalid normalization method. Supported methods are 'RLE' and 'TMM'")
  }
  
  # Combine gene_info and normalizedCount
  data2 <- cbind(gene_info, normalizedCount)
  
  write.table("#1.2", filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(paste(nrow(data2), ncol(data2) - 2, sep = "\t"), filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(data2, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE)
}
