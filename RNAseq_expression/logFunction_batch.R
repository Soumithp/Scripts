args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3 ) {
  stop("Usage: Rscript --vanilla logFunction_batch.R <script> <inputFile> <logBase>", call.=FALSE)
} else if (length(args) == 3) {
  script <- args[1]
  inputFile <- args[2]
  logBase <- as.numeric(args[3])
}

source(script)

logGct(inputFile, logBase)
