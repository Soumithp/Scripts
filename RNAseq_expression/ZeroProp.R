####auhtor: soumith paritala
####                    modifications                                                 date:
#         added intermediate file for cross checking the result                   jan-15-2024
#         removed ceiling on the accepted number & added "=" to zerothresh        feb-13-2024





ZeroProp <- function(
    inputFile = "",
    zeroProp = 0.9
) {
  data <- read.delim(inputFile, header = TRUE, skip = 2)
  NoCols <- ncol(data) - 2
  zeroProp <- as.numeric(zeroProp)
  acceptedNum <- NoCols * zeroProp
  zeroCounts <- rowSums(data == 0)
  intermediateResult <- cbind(data[, c(1, 2)], ZeroCount = zeroCounts)
  
  fileName <- tools::file_path_sans_ext(inputFile)
  outputFile <- paste0(fileName, "_0p_total.txt")
  write.table(intermediateResult, outputFile, quote = FALSE, sep = "\t", row.names = FALSE)
  
  zeroThresh <- data[rowSums(data == 0) <= acceptedNum, ]
  zeroPercent<-zeroProp*100
  
  fileName <- tools::file_path_sans_ext(inputFile)

  geneInfo <- data[, c(1, 2)]
  zeroThresh <- merge(geneInfo, zeroThresh, by = c(1, 2))
  
  # # Add a new column for values of genes in the subset
  # intermediateResult$ValuesAfterCutoff <- NA
  # idx <- match(intermediateResult[, c(1, 2)], zeroThresh[, c(1, 2)])
  # intermediateResult$ValuesAfterCutoff[idx] <- zeroThresh$YourColumnName
  

  FileName <- paste0(fileName,"_0p",zeroPercent,".gct",sep="")
  
  write.table( "#1.2", FileName, quote=F,sep="\t",row.names=F,col.names=F)
  write.table( paste( nrow(zeroThresh), ncol(zeroThresh)-2, sep="\t"), FileName, quote=F, sep="\t", row.names=F, col.names=F, append=T)
  write.table( zeroThresh, FileName, quote=F, sep="\t", row.names=F, col.names=T, append=T)
  
  return(0)
}
  
  