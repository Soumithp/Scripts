VarFilter <- function(
    inputFile,
    var= "cv",
    var.cutoff = 0.1
) {
  fileName <- tools::file_path_sans_ext(inputFile)
  inputFile <- read.delim(inputFile, header = TRUE, skip = 2)
  row_std <- apply(inputFile[c(-1, -2)], 1, sd)
  row_mean <- apply(inputFile[c(-1, -2)], 1, mean)
  MeanAndStD <- cbind(inputFile, mean = row_mean, std = row_std)
  MeanStDCV <- cbind(MeanAndStD, cv = MeanAndStD$std / MeanAndStD$mean)

  IntermediateFileName <- paste(fileName, "_mean_std_CV.txt", sep = "")
  write.table(MeanStDCV[, c(1, 2, ncol(MeanStDCV) - 2, ncol(MeanStDCV) - 1, ncol(MeanStDCV))], 
              IntermediateFileName, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    if (var == "cv") {
    FinalOutput <- subset(MeanStDCV, cv > var.cutoff)
  } else if (var == "std") {
    FinalOutput <- subset(MeanStDCV, std < var.cutoff)
  } else {
    stop("Invalid value for 'var'. Use 'cv' or 'std'.")
  }
  
  FinalOutput <- FinalOutput[, !(names(FinalOutput) %in% c("cv", "std", "mean"))]
  
  FileName <- paste(fileName, "_", var, var.cutoff, ".gct", sep="")
  
  write.table("#1.2", FileName, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(paste(nrow(FinalOutput), ncol(FinalOutput) - 2, sep = "\t"), FileName, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(FinalOutput, FileName, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE)
  
  return(0)

  
}
