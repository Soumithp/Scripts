# IntersampleCorrelation
#                            09/17/2015
#        update output.name  10/07/2015
#        input.exp.filename removed from output.name 03/15/2016
#        +1 option added  01/31/2023
#
#  input: gene expression data (.gct)
#   geo.mean is calclutaed on log

IntersampleCorrelation <- function(
    input.exp.filename,
    reference.method="median", # or "mean", "geo.mean"
    add.one="T", # +1 for log
    output.name="IntCor"
    )
{
    d <- read.delim(input.exp.filename,header=T,skip=2,check.names=F)
    d <- d[-c(1:2)]
    sample.name <- colnames(d)

    if (add.one=="T"){
        d1 <- d+1
    }else{
        d1 <- d
    }

    if (reference.method=="median"){
        ref <- apply(d1,1,median,na.rm=T)
    }
    if (reference.method=="mean"){
        ref <- apply(d1,1,mean,na.rm=T)
    }
    if (reference.method=="geo.mean"){
        ref <- apply(log(d1),1,prod,na.rm=T)^(1/dim(d)[2])
    }

    pearson.r <- cor(d,ref,use="complete.obs",method="pearson")
    spearman.rho <- cor(d,ref,use="complete.obs",method="spearman")

    pearson.r.on.log <- cor(log(d1),log(ref),use="complete.obs",method="pearson")
    spearman.rho.on.log <- cor(log(d1),log(ref),use="complete.obs",method="spearman")

    output <- cbind(sample.name,pearson.r,spearman.rho,pearson.r.on.log,spearman.rho.on.log)
    colnames(output) <- c("sample.name","pearson.r","spearman.rho","pearson.r.on.log","spearman.rho.on.log")
    output.name <- paste(output.name,"_",reference.method,".txt",sep="")
#    output.name <- paste(output.name,"_",gsub(".gct","",input.exp.filename),"_",reference.method,".txt",sep="")

    write.table(output,output.name,quote=F,sep="\t",row.names=F)

}




