##Normalization batch code after generating rawcounts
##date_modified:19-jan-2024
##updated: made a single file for normalization
##input = rawcounts .gct file , normalization method
##output= normalized .gct file 
##input file, source code and batch code should be in the same directory

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript --vanilla normalizeData_batch.R <script> <inputFile> <normalizationMethod>", call.=FALSE)
} else if (length(args) == 3) {
  script <- args[1]
  inputFile <- args[2]
  normalizationMethod <- args[3]
  
}

# current_dir <- getwd()
# setwd(current_dir)
source(script)


# inputFile <- ".gct"
# normalizationMethod <- "RLE"

normalizeData(inputFile, normalizationMethod)
