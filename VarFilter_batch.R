args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4 ) {
  stop("Usage: Rscript --vanilla VarFilter_batch.R <script> <inputFile> <var> <varCutoff>", call.=FALSE)
} else if (length(args) == 4) {
  script <- args[1]
  inputFile <- args[2]
  var <- args[3]
  varCutoff <- as.numeric(args[4])
}

# current_dir <- getwd()
# setwd(current_dir)
source(script)

# inputFile <- "DENCCI4_7rh_RLE_symbolmed_humanised_Tumor.gct"
# var = "cv"  # Choose "cv" or "std"
# varCutoff = 0.5  # Set your desired cutoff value

VarFilter(inputFile, var, varCutoff)
