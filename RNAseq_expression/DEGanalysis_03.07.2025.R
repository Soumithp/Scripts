
setwd("C:/Users/s226953/OneDrive - University of Texas Southwestern/Desktop/PLScp_EGCG_CVC_05.30.2024/")
##############################################################################
# 1) Load libraries
##############################################################################
# install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)   # optional for plotting, if desired

read_gct <- function(gct_path) {
  con <- file(gct_path, "r")
  gct_version <- readLines(con, n=1)  # e.g. "#1.2"
  dim_line <- readLines(con, n=1)
  
  dims   <- strsplit(dim_line, "\t")[[1]]
  n_rows <- as.integer(dims[1])
  n_cols <- as.integer(dims[2])
  
  all_lines <- readLines(con)
  close(con)
  
  tf <- textConnection(all_lines, "r")
  df <- read.table(tf,
                   header=TRUE,
                   sep="\t",
                   quote="",
                   comment.char="",
                   check.names=FALSE,
                   nrows=n_rows)
  close(tf)
  
  # Typically: first two columns are "Name" and "Description"
  colnames(df)[1:2] <- c("Name", "Description")
  
  sample_cols <- colnames(df)[3:(2 + n_cols)]
  count_mat <- as.matrix(df[, sample_cols, drop=FALSE])
  rownames(count_mat) <- df$Name
  
  # Remove all-zero rows
  keep <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep, , drop=FALSE]
  
  list(counts = count_mat, samples = sample_cols)
}

##############################################################################
# 3.1) Read GCT: raw counts
##############################################################################
gct_data <- read_gct("C:/Users/s226953/OneDrive - University of Texas Southwestern/Desktop/PLScp_EGCG_CVC_05.30.2024/PLScp_EGCG_CVC_Rawcountsmatrix.gct")
counts   <- gct_data$counts  # Genes x Samples

##############################################################################
# 3.2) Read metadata
#     Expect columns: Sample, sex, treatment
##############################################################################
meta_data <- read.table("C:/Users/s226953/OneDrive - University of Texas Southwestern/Desktop/PLScp_EGCG_CVC_05.30.2024/metadata.txt",
                        header=TRUE,
                        sep="\t",
                        stringsAsFactors=FALSE)

# Make rownames = Sample
rownames(meta_data) <- meta_data$Sample

# Check sample names match
if (!all(colnames(counts) %in% rownames(meta_data))) {
  stop("Some GCT columns not found in metadata's Sample column.")
}

##############################################################################
# 3.3) Function to run DESeq2 for one sex
#     - Subset counts & metadata by that sex
#     - Relevel baseline to "CDA"
#     - Compare each other group vs. CDA
#     - Output results & plot MA
##############################################################################
run_deseq_for_sex <- function(sex_label, counts, meta_data, out_prefix="DEG_results") {
  # Subset metadata for this sex
  meta_sub <- meta_data[meta_data$sex == sex_label, ]
  
  # Subset counts
  sel_samples <- meta_sub$Sample
  counts_sub  <- counts[, sel_samples, drop=FALSE]
  
  # Build DESeq dataset
  dds <- DESeqDataSetFromMatrix(countData=counts_sub,
                                colData=meta_sub,
                                design=~ treatment)
  
  # Relevel so "CDA" is baseline
  # (If your baseline is different, change it here)
  dds$treatment <- factor(dds$treatment)
  dds$treatment <- relevel(dds$treatment, ref="CH")
  
  # Run DESeq
  dds <- DESeq(dds)
  
  # Identify all treatment levels
  trt_levels <- levels(dds$treatment)
  # We'll do each level vs "CDA"
  baseline   <- "CH"
  # The list of comparisons is all levels except "CDA"
  compare_levels <- setdiff(trt_levels, baseline)
  
  # For each comparison, write results & do an MA plot
  for (grp in compare_levels) {
    # e.g., CH vs CDA
    res <- results(dds, contrast=c("treatment", grp, baseline))
    
    # Order by padj
    res <- res[order(res$padj), ]
    
    # Summaries
    cat("\n===== ", sex_label, ": ", grp, " vs ", baseline, " =====\n")
    summary(res)
    cat("padj < 0.1: ", sum(res$padj < 0.1, na.rm=TRUE), "\n")
    
    # Plot MA
    plotMA(res, main=paste0("DESeq2 ", sex_label, ": ", grp, " vs ", baseline),
           ylim=c(-2,2))
    
    # Filter padj < 0.1
    res_filt <- subset(res, padj < 0.1)
    
    # Save table
    out_file <- paste0(out_prefix, "_", sex_label, "_vs", baseline, "_", grp, "_padj0.1.txt")
    write.table(as.data.frame(res_filt), file=out_file,
                row.names=TRUE, quote=FALSE, sep="\t")
    
    cat("Wrote ", out_file, " with ", nrow(res_filt), " genes.\n")
  }
  
  return(dds)
}

##############################################################################
# 3.4) Run for Males & Females
##############################################################################
dds_male <- run_deseq_for_sex(
  sex_label="M",
  counts=counts,
  meta_data=meta_data,
  out_prefix="DEG_results_male"
)

dds_female <- run_deseq_for_sex(
  sex_label="F",
  counts=counts,
  meta_data=meta_data,
  out_prefix="DEG_results_female"
)


