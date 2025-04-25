setwd("C:/Users/s226953/OneDrive - University of Texas Southwestern/Desktop/PLScp_EGCG_CVC_05.30.2024/")

# Load libraries
library(DESeq2)
library(data.table)

source("DEG_analysis_function.R")

# Example usage:
gct_file <- "C:/Users/s226953/OneDrive - University of Texas Southwestern/Desktop/PLScp_EGCG_CVC_05.30.2024/PLScp_EGCG_CVC_Rawcountsmatrix.gct"  
metadata_file <- "C:/Users/s226953/OneDrive - University of Texas Southwestern/Desktop/PLScp_EGCG_CVC_05.30.2024/metadata.txt"


# Run the analysis
dds <- runDEGanalysis(
  gct_file = gct_file,
  metadata_file = metadata_file,
  output_prefix = "deg_results"
)

