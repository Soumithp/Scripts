wkdir = "/work/SCCC/s184942/LOOCV.Survival/CGMH_validation/Roessler_cirrhosis_wo.w"
setwd(wkdir)
source("NTPtoKMsummary.R")
library(foreach)
library(doParallel)

cl <-detectCores()
registerDoParallel(cl) 


setwd(paste0(wkdir, "/signatures"))
lst.sig <- list.files(pattern=".txt")
setwd(wkdir)
dir.create("output")
dir.create("output/rec")
dir.create("output/death")
dir.create("output/late.rec")
foreach (i = seq_along(lst.sig))%dopar%{
  cond <- gsub("_LOOCVSurvivalSignatures.txt", "", gsub("Variation.Filter_", "", lst.sig[i]))
  sig <- read.delim(paste0("signatures/",lst.sig[i] ),header =T)
  if(length(unique(sig[,3]))==1){
    print(paste0(cond, " has no proper signature"))
  }else{
  NTPtoKMsummary(
  expression.data = "ROESSLER_HEPATOCELLULAR_CARCINOMA_CHINA_LIVER_symbolmed_cirrhosis.gct",
  input.clinical.data = "ROESSLER_HEPATOCELLULAR_CARCINOMA_CHINA_LIVER_outcome_cirrhosis.txt",
  signature.table = paste0("signatures/", lst.sig[i]),
  pt.ID = "sample.names", ##please set the patient IDs which are used in gct
  ntp.num.resamplings = 1000, #NTP permutation
  ntp.output.name = paste0("output/rec/RoesslerLiver.rec_CMGHsig_", cond), 
  FDR.sample.bar=0.25,
  
  weight.genes="F",
  within.sig = "F",
  ##KM setting
  time ="time.to.rec",  #input the column name of time field
  censor = "rec",  #input the column name of censor field
  conf.int = FALSE,
  mark.time = FALSE,
  fun = "event",
  ylim = c(0,1)
  
)
  }




  cond <- gsub("_LOOCVSurvivalSignatures.txt", "", gsub("Variation.Filter_", "", lst.sig[i]))
  sig <- read.delim(paste0("signatures/",lst.sig[i] ),header =T)
  if(length(unique(sig[,3]))==1){
    print(paste0(cond, " has no proper signature"))
  }else{
    
   NTPtoKMsummary(
    expression.data = "ROESSLER_HEPATOCELLULAR_CARCINOMA_CHINA_LIVER_symbolmed_cirrhosis.gct",
    input.clinical.data = "ROESSLER_HEPATOCELLULAR_CARCINOMA_CHINA_LIVER_outcome_cirrhosis.txt",
    signature.table = paste0("signatures/", lst.sig[i]),
    pt.ID = "sample.names", ##please set the patient IDs which are used in gct
    ntp.num.resamplings = 1000, #NTP permutation
    ntp.output.name = paste0("output/death/RoesslerLiver.death_CMGHsig_", cond), 
    FDR.sample.bar=0.25,
    
    weight.genes="F",
    within.sig = "F",
    ##KM setting
    time ="time.to.last",  #input the column name of time field
    censor = "death",  #input the column name of censor field
    conf.int = FALSE,
    mark.time = FALSE,
    fun = "event",
    ylim = c(0,1)
    
  )
  }





  cond <- gsub("_LOOCVSurvivalSignatures.txt", "", gsub("Variation.Filter_", "", lst.sig[i]))
  sig <- read.delim(paste0("signatures/",lst.sig[i] ),header =T)
  if(length(unique(sig[,3]))==1){
    print(paste0(cond, " has no proper signature"))
  }else{
    
    NTPtoKMsummary(
      expression.data = "ROESSLER_HEPATOCELLULAR_CARCINOMA_CHINA_LIVER_symbolmed_cirrhosis.gct",
      input.clinical.data = "ROESSLER_HEPATOCELLULAR_CARCINOMA_CHINA_LIVER_outcome_cirrhosis.txt",
      signature.table = paste0("signatures/", lst.sig[i]),
      pt.ID = "sample.names", ##please set the patient IDs which are used in gct
      ntp.num.resamplings = 1000, #NTP permutation
      ntp.output.name = paste0("output/late.rec/RoesslerLiver.late.rec_CMGHsig_", cond), 
      FDR.sample.bar=0.25,
      
      weight.genes="F",
      within.sig = "F",
      ##KM setting
      time ="time.to.rec",  #input the column name of time field
      censor = "late.rec",  #input the column name of censor field
      conf.int = FALSE,
      mark.time = FALSE,
      fun = "event",
      ylim = c(0,1)
      
    )
  }
}


stopCluster()
