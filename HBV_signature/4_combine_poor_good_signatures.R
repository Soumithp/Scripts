#############################################################
# combine_poor_good_signatures.R      Jul.25,2012
#                        (Kensuke Kojima)
#
# input
#(1)output file of survivalgene.R: .txt file
#(2)output file of loocv.R (good.features): .txt file
#(3)output file of loocv.R (poor.features): .txt file
#
# output: gene.signatures.w/statistic.values: .txt
#
# parameters
#(1) cutoff.percent:
#   Signature genes chosen in the leave-one-out trials equal or greater than this cut-off (in %)
#############################################################

combine_poor_good_signatures <- function(
    input.file.surv,
    input.file.poor,
    input.file.good,
    cutoff.percent = 100,
    output.file = "LOOCV_Survival_Signatures"
    )
{

#define.data.type
 cutoff.percent <-  as.numeric(cutoff.percent)
 cutoff.percent <- cutoff.percent/100

#imput.files
 surv <- read.delim(input.file.surv,header=T)
 loocv.poor <- read.delim(input.file.poor,header=T)
 loocv.good <- read.delim(input.file.good,header=T)

 num.genes <- dim(surv)[1]
 num.samples <- dim(loocv.poor)[2]-4


### poor.features
 poor.features <- merge(surv,loocv.poor,all=F)
 poor.features.essense <- poor.features[c(1:5,num.samples+9)]

 #insert.flag(direction)
  annot <- poor.features.essense[c(1,2)]
  poor.features.essense <- poor.features.essense[-c(2)]

  poor1_good2 <- matrix(1,num.genes,1)
  annot.poor1_good2<- cbind(annot,poor1_good2)
  poor.features.essense <- merge(annot.poor1_good2,poor.features.essense)

 #filtering
  poor.features.essense.filtered <- poor.features.essense[which(poor.features.essense$percent.poor>=cutoff.percent),]

 #change.column.name
  names(poor.features.essense.filtered)[7] <- "percent"



### good.features
 good.features <- merge(surv,loocv.good,all=F)
 good.features.essense <- good.features[c(1:5,num.samples+9)]

  #insert.flag(direction)
  annot <- good.features.essense[c(1,2)]
  good.features.essense <- good.features.essense[-c(2)]

  poor1_good2 <- matrix(2,num.genes,1)
  annot.poor1_good2 <- cbind(annot,poor1_good2)
  good.features.essense <- merge(annot.poor1_good2,good.features.essense)

 #filtering
  good.features.essense.filtered <- good.features.essense[which(good.features.essense$percent.good>=cutoff.percent),]
 #change.column.name
  names(good.features.essense.filtered)[7] <- "percent"




### merge_poor.good.features
 poor.good.merged.list <- rbind(poor.features.essense.filtered,good.features.essense.filtered)

 sortlist3 <- order(poor.good.merged.list$statistic,decreasing=T)
 final.list <- poor.good.merged.list[sortlist3,]
 final.list <- final.list[-c(5:7)]

### output.file
 output.file <- paste(output.file,".txt",sep="")
 write.table(final.list,output.file,quote=F,sep="\t",row.names=F,col.names=c("identifier","description","poor1_good2","cox.score"))
}
