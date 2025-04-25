args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2 ) {
  stop("Usage: Rscript --vanilla normalize_mutationalpatterns.R <in> <out>", call.=FALSE)
} else if (length(args) == 2) {
  in_file <- args[1]
  out_file <- args[2]
}

in_mat <- read.csv(in_file, row.names = 1)
norm_mat <- t(t(in_mat)/colSums(in_mat))
write.csv(norm_mat, file = out_file, row.names = TRUE)