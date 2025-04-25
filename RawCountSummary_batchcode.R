##Summarization step of the raw counts from feature counts output
##step:3(summarizing the counts of all samples) from PipelineSOP01_RNAseq_expression_2023.12.14
##SOP version: 1.0

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript --vanilla RawCountSummary_batchcode.R <>", call.=FALSE)
} else if (length(args) >= 4) {
  script <- args[1]
  projectName <- args[2]
  annot <- args[3]
  featureCountsFiles <- args[4:length(args)]
  
}

# outputDir="path/to/output"
# projectName = "Sample_project_Name" 
gtfAnnot = read.table(annot, sep="\t", header=T)

source(script)

#rawcounts_summary
#provide the path to the featurecounts.primary.txt files 
RawCountSummary(files = featureCountsFiles, 
  outputFileName = projectName, 
  gtfAnnot = gtfAnnot)
# rawcountsummary(cohort_list_file= "path/to/the/cohort_list", directory = "/path/to/the/folder/", outputFileName = projectName, gtfAnnot=gtfAnnot)
