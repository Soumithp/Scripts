rm(list=ls())

###Parameters##
root.dir <- '/work/SCCC/s184942/LOOCV.Survival/CGMH_log';setwd(root.dir)
expression.data <- "HBV_CGMH_renormalized_symbolmed_protein.coding_log2.gct"
summary.table <- "HBV_CGMH_renormalized_symbolmed_log2_protein.coding_summary.txt"
cv.c <- c(0.1,0.5,1,2) # cv you want to try
diff.c <- c(0.5,1,2)
non0per <- c(0.1,0.2,0.5,0.9)
input.clinical.data = 'CGMH_HBVcir_n89_2019.12.12.txt'
pt.ID <- "sample.names" ##please set the patient IDs which are used in gct
nperm <- 1000 ##SurvivalGene permutation
nresample <- 1000 ##LOOCV number of resampling

ntp.num.resamplings <- 1000 #NTP permutation
##KM setting
time ="days.to.hcc"  #input the column name of time field
censor = "hcc"
conf.int <- FALSE
mark.time <- FALSE
fun <-"event"

################


#CV
library(foreach)
library(doParallel)
num.cores <- detectCores() - 1 
cl <- makeCluster(num.cores)
registerDoParallel(cl)
dir.create("signatures")
dir.create("SurvivalGene")

for( i in seq_along(cv.c)){
  for (j in seq_along(diff.c)){
    foreach(k = seq_along(non0per)) %dopar% {
  

  source("2_SurvivalGene.R")
  source("3_loocv.nn_cox_GP.R")
  source("4_combine_poor_good_signatures.R")
  source("NTP.R")  
  
dir.create(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]))


exp <- read.delim(expression.data,skip=2, header =T)
smr <- read.delim(summary.table, header =T)

smr1 <- subset(smr, cv>cv.c[i]&diff>diff.c[j]&num.zero.per<non0per[k])
selected.genes <- smr1[,1]
exp1 <- exp[which(exp[,1]%in%selected.genes),]
gene.num <- nrow(exp1)
sample.num <- ncol(exp1)-2

if(gene.num>50){
write.table("#1.2", paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "gct",sep = "."),sep = "\t", col.names = F, row.names = F, quote =F)
write.table(data.frame(gene.num, sample.num), paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "gct",sep = "."),sep = "\t", col.names = F, row.names = F, quote =F,append=T)
write.table(exp1, paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "gct",sep = "."),sep = "\t", col.names = T, row.names = F, quote =F,append=T)


file.copy(from=paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "gct",sep = "."), 
          to=paste(gsub("'", "", root.dir), paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "/"), 
          copy.mode = TRUE)
file.remove(paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "gct",sep = "."))

###SurvivalGene 

SurvivalGene(
  paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "gct",sep = "."),sep = "/"),
  input.filename.clinical = input.clinical.data, 
  output.file=paste(paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), sep = "/"), "SurvivalGene",sep = "_"), 
  
  time.field=time,##PLEASE CHANGE
  censor.field=censor,##PLEASE CHANGE
  
  statistic.selection="cox.score",
  trim.percent.2.side=0,
  nperm=nperm,
  rnd.seed=1234567,
  emp.stat.dist="T"
)

file.copy(from=paste(paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), sep = "/"), "SurvivalGene.txt",sep = "_"), 
          to=paste(gsub("'", "", root.dir), "SurvivalGene", sep = "/"), 
          copy.mode = TRUE)


##LOOCV


loocv.nn_cox_GP(
  paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "gct",sep = "."),sep = "/"),
  input.clinical.data,
  
  # for CoxScore
   time.field=time,##PLEASE CHANGE
  censor.field=censor,##PLEASE CHANGE
  
  trim.percent.2.side=0,
  cox.score.sig=0.05,  # 2-sided p-value. >=1 is regarded as # of marker genes
  # assume normal dist
  emp.cox.file=paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "SurvivalGene_emp.stat.txt",sep = "_"),sep = "/"),   # emripical dist of Cox score to standardize
  
  # Prediction
  temp.nn.wt="T",
  dist.selection="cosine",  # "correlation" or "cosine"
  
  # number of resampling to generate null dist for prediction stat
  nresample=nresample,
  
  # random seed
  rnd.seed=4675921,
  
  output.file=paste(paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), sep = "/"), "LOOCVSurvival",sep = "_")
)


##combine_poor_good_signatures


  combine_poor_good_signatures (
    paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "SurvivalGene.txt",sep = "_"),sep = "/"),
    paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "LOOCVSurvival_poor.features.txt",sep = "_"),sep = "/"),
    paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "LOOCVSurvival_good.features.txt",sep = "_"),sep = "/"),
    output.file = paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste(paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), "LOOCVSurvivalSignatures",sep = "_"),sep = "/")
    )
  
  
  sig <- read.delim(paste(paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), sep = "/"), "LOOCVSurvivalSignatures.txt",sep = "_"),header =T)
  
  if(length(unique(sig[,3]))==2){
  file.copy(from=paste(paste(paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), paste("Variation.Filter", paste0("cv", cv.c[i],"_diff", diff.c[j], "_no0per", non0per[k]), sep = "_"), sep = "/"), "LOOCVSurvivalSignatures.txt",sep = "_"),
            to = "signatures", 
            copy.mode = TRUE)
                                }
  
    }
  }
 }
}

stopCluster(cl)


