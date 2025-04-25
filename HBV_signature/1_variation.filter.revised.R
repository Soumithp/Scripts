#############################
# Valiation.Filter.R
#  Xintong Chen  03.19.2013
#
#############################
Variation.Filter <- function(
    input.exp.filename="",
    ordering.method="CV",#or "STD" or "MAD"
    cutoff.method="probe.number",# "probe.number" or "cutoff.value"
    cutoff.value=5000,
    output.filename="Variation.Filter"
    )
{


#for GenePattern
cutoff.value <- as.numeric(cutoff.value)

#input.exp.file
d <- read.delim(input.exp.filename,header=T,skip=2)
d.data <- d[-c(1,2)]
d.probe <- d[c(1,2)]

#calcurate s.d. and mean for each probe.
std <- apply(d.data,1,sd,na.rm=T)
m <- apply(d.data,1,mean,na.rm=T)
mads<- apply(d.data,1,mad,na.rm=T)

if(ordering.method=="CV")
{
cv <- std/m
output.data.cv <- cbind(d.probe,d.data,cv)
write.table (output.data.cv,paste(output.filename,"_cv.txt",sep=""),quote=F,sep="\t",row.names=F)
sort.cv <- output.data.cv[order(output.data.cv$cv,decreasing=T),]

  if(cutoff.method == "cutoff.value"){
      sort.filtered.cv <- subset(sort.cv,cv>cutoff.value)
      num.col.cv <- dim(sort.filtered.cv)[2]
      sort.filtered.cv <- sort.filtered.cv[-c(num.col.cv)]
      num.col.cv <- dim(d.data)[2]
      num.row.cv <- dim(sort.filtered.cv)[1]
      write.table("#1.2",paste(output.filename,"_cv",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
      write.table(paste(num.row.cv,num.col.cv,sep="\t"),paste(output.filename,"_cv",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
      write.table(sort.filtered.cv,paste(output.filename,"_cv",cutoff.value,".gct",sep=""),
              quote=F,sep="\t",row.names=F,col.names=T,append=T)
  }

  if(cutoff.method == "probe.number"){
       if(cutoff.value<1){
          stop("###cutoff.value must be integer!###")
      }
      sort.filtered.cv <- head(sort.cv,n=cutoff.value)
      num.col.cv <- dim(sort.filtered.cv)[2]
      sort.filtered.cv <- sort.filtered.cv[-c(num.col.cv)]
      num.col.cv <- dim(sort.filtered.cv)[2]-2
      num.row.cv <- dim(sort.filtered.cv)[1]
      write.table("#1.2",paste(output.filename,"_cv",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
      write.table(paste(num.row.cv,num.col.cv,sep="\t"),paste(output.filename,"_cv",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
      write.table(sort.filtered.cv,paste(output.filename,"_cv",cutoff.value,".gct",sep=""),
              quote=F,sep="\t",row.names=F,col.names=T,append=T)
  }
}

if(ordering.method=="STD")
{
output.data.std <- cbind(d.probe,d.data,std)
write.table (output.data.std,paste(output.filename,"_sd.txt",sep=""),quote=F,sep="\t",row.names=F)
sort.std <- output.data.std[order(output.data.std$std,decreasing=T),]

 if(cutoff.method == "cutoff.value"){
      sort.filtered.std <- subset(sort.std,std>=cutoff.value)
      num.col.std <- dim(sort.filtered.std)[2]
      sort.filtered.std <- sort.filtered.std[-c(num.col.std)]
      num.col.std <- dim(sort.filtered.std)[2]-2
      num.row.std <- dim(sort.filtered.std)[1]
      write.table("#1.2",paste(output.filename,"_sd",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
      write.table(paste(num.row.std,num.col.std,sep="\t"),paste(output.filename,"_sd",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
      write.table(sort.filtered.std,paste(output.filename,"_sd",cutoff.value,".gct",sep=""),
              quote=F,sep="\t",row.names=F,col.names=T,append=T)
  }
 if(cutoff.method == "probe.number"){
      sort.filtered.std <- head(sort.std,n=cutoff.value)
      num.col.std <- dim(sort.filtered.std)[2]
      sort.filtered.std <- sort.filtered.std[-c(num.col.std)]
      num.col.std <- dim(sort.filtered.std)[2]-2
      num.row.std <- dim(sort.filtered.std)[1]
      write.table("#1.2",paste(output.filename,"_sd_",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
      write.table(paste(num.row.std,num.col.std,sep="\t"),paste(output.filename,"_sd",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
      write.table(sort.filtered.std,paste(output.filename,"_sd",cutoff.value,".gct",sep=""),
              quote=F,sep="\t",row.names=F,col.names=T,append=T)
  }
}




if(ordering.method=="MAD")
{
  output.data.mad <- cbind(d.probe,d.data,mads)
  write.table (output.data.mad,paste(output.filename,"_MAD.txt",sep=""),quote=F,sep="\t",row.names=F)
  sort.mad <- output.data.mad[order(output.data.mad$mads,decreasing=T),]
  
  if(cutoff.method == "cutoff.value"){
    sort.filtered.mad <- subset(sort.mad,mads>=cutoff.value)
    num.col.mad <- dim(sort.filtered.mad)[2]
    sort.filtered.mad <- sort.filtered.mad[-c(num.col.mad)]
    num.col.mad <- dim(sort.filtered.mad)[2]-2
    num.row.mad <- dim(sort.filtered.mad)[1]
    write.table("#1.2",paste(output.filename,"_mad",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    write.table(paste(num.row.mad,num.col.mad,sep="\t"),paste(output.filename,"_mad",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table(sort.filtered.mad,paste(output.filename,"_mad",cutoff.value,".gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=T,append=T)
  }
  if(cutoff.method == "probe.number"){
    sort.filtered.mad <- head(sort.mad,n=cutoff.value)
    num.col.mad <- dim(sort.filtered.mad)[2]
    sort.filtered.mad <- sort.filtered.mad[-c(num.col.mad)]
    num.col.mad <- dim(sort.filtered.mad)[2]-2
    num.row.mad <- dim(sort.filtered.mad)[1]
    write.table("#1.2",paste(output.filename,"_mad_",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    write.table(paste(num.row.mad,num.col.mad,sep="\t"),paste(output.filename,"_mad",cutoff.value,".gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table(sort.filtered.mad,paste(output.filename,"_mad",cutoff.value,".gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=T,append=T)
  }
}
}

