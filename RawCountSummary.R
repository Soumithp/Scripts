##Summarization step
##step:3(summarizing the counts of all samples) from PipelineSOP01_RNAseq_expression_2024.01.08.docx
### modifications                               date:
# changed to get only gene_id & gene_name       Feb-9th-2024



RawCountSummary <- function(files, outputFileName, gtfAnnot) {
  # setwd(directory)
  # samples <- dir()
  # print(samples) ## samples[! names(samples) %in% c('multiqc_report.html','multiqc_data' or from the cohort list)]
  # gtfAnnot$gene_id = sub("\\..*$","",gtfAnnot$gene_id)

  samples <- sub(".featureCounts.primary.txt", "", basename(files))

  count1 <- lapply(1:length(samples), function(i) {
    dd <- samples[i]
    file <- files[i]
    cat(dd,"\n")
    read.table(file, sep = "\t", header = TRUE)
  })
  
  # countMatrix <- do.call(cbind, lapply(count1, function(x) x[, 7]))
  # gene_id <- sapply(strsplit(as.character(count1[[1]]$Geneid), "[.]"), function(x) x[1])
  gene_id <- as.character(count1[[1]]$Geneid)
  count2 <- do.call(cbind, lapply(count1, function(x) x[, 7]))
  colnames(count2) <- samples
  rownames(count2) <- gene_id
  count3 = data.frame(gene_id,count2)
  count = merge(gtfAnnot,count3,by="gene_id")
  
 
  count2 <- count[, !(names(count) %in% c("chr", "strand", "start", "end", "gene_type"))]
  
  
  
  filename= paste0(outputFileName,"_Rawcountsmatrix.gct")
                  
  write.table( "#1.2", filename, quote=F,sep="\t",row.names=F,col.names=F)
  write.table( paste( nrow(count2), ncol(count2)-2, sep="\t"), filename, quote=F, sep="\t", row.names=F, col.names=F, append=T)
  write.table( count2, filename, quote=F, sep="\t", row.names=F, col.names=T, append=T)
 
}

