#setwd("/work/SCCC/s184942/LOOCV.Survival/CMGH_log")
summary_exp=function(exp){
  #input.exp.file
  d <- exp#read.delim("CGMH_HBV_symbolmed_cv0.01_no0in90per.gct",header=T,skip=2)
  d.data <- d[-c(1,2)]
  d.probe <- d[c(1,2)]
  
  #d.data <- log2(d.data+1)
  
  
  d.probe$std <- apply(d.data,1,sd,na.rm=T)
  d.probe$m <- apply(d.data,1,mean,na.rm=T)
  d.probe$median <- apply(d.data,1,median,na.rm=T)
  
  d.probe$mads<- apply(d.data,1,mad,na.rm=T)
  d.probe$max <- apply(d.data,1, max, na.rm=T)
  d.probe$min <- apply(d.data,1, min, na.rm=T)
  d.probe$cv <- d.probe$std/d.probe$m
  d.probe$cv.abs <- abs(d.probe$cv)
  d.probe$diff <- d.probe$max-d.probe$min
  d.probe$num.zero.per <- rowSums(d.data==0)/ncol(d.data)
  
  #write.table(d.probe, "CGMHn82_RLE_symbolmed_cv0.01_no0in90per_summary.txt",sep ="\t",col.names = T, row.names = F,quote =F)
  return(d.probe)
}


