args = commandArgs(trailingOnly=TRUE)
if (length(args) !=3 ) {
  stop("Usage: Rscript --vanilla ZeroProp_batchcode.R <script> <inputFile> <zeroProp>", call.=FALSE)
} else if (length(args) == 3) {
  script <- args[1]
  inputFile <- args[2]
  zeroProp <- args[3]
  
}

source(script)

# inputFile <- "DENCCI4_7rh_RLE_symbolmed_humanised_Tumor.gct"
# zeroProp <- 0.9

ZeroProp(inputFile, zeroProp)
