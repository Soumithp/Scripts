logGct <- function(inputFile, logBase = 2, add.num= 1) {
  
  normalizedCount <- read.delim(inputFile, header = TRUE, skip = 2)
  nonCountCols <- c(1:2)
  fileName <- tools::file_path_sans_ext(inputFile)
  nonCountColsData <- normalizedCount[, nonCountCols, drop = FALSE]
  otherColsData <- normalizedCount[, -nonCountCols, drop = FALSE]

  logTransformedData <- log(otherColsData + add.num, base = logBase)

  FinalOutput <- cbind(nonCountColsData, logTransformedData)

  FileName <- paste0(fileName, paste0("_log", logBase), ".gct")
  
  # Write the result to a file
  write.table("#1.2", FileName, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(paste(nrow(FinalOutput), ncol(FinalOutput) - 2, sep = "\t"), FileName, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(FinalOutput, FileName, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE)
}
