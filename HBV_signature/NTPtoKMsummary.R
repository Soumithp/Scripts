######################################
#
#   ntp_toKMsummary
#   5/9/2018, summary table update


NTPtoKMsummary <- function(
expression.data = "LuminexPLS_Lgh001_1to3Candidate_serum_n79_FI_normport_LR.gct",
input.clinical.data = "LGh001_clinical_outcome_n79_forSerum.txt",
signature.table = "Luminex_PLS4cohort_es.ex_w.w_forNTP.txt",
pt.ID = "bkg_ID", ##please set the patient IDs which are used in gct
ntp.num.resamplings = 1000, #NTP permutation
ntp.output.name = "Luminex_PLS_w.w_Lgh001n79_rec", 
FDR.sample.bar=0.25,

weight.genes="T",
within.sig="F", 
##KM setting
time ="days.to.rec",  #input the column name of time field
censor = "rec",  #input the column name of censor field
conf.int = FALSE,
mark.time = FALSE,
fun = "event", #for KM
ylim = c(0,1)
){
#################################################
###############
#     NTP     #
###############

  try(NTP(
    expression.data,
    signature.table,
    output.name = ntp.output.name,
    
    distance.selection="cosine", # "correlation" or "cosine"
    weight.genes=weight.genes,   # only for 2 cls
    
    num.resamplings=ntp.num.resamplings,
    within.sig=within.sig,
    
    GenePattern.output="T",
    signature.heatmap="T",
    p.sample.bar=0.05, # NA if not needed
    FDR.sample.bar=FDR.sample.bar, # NA if not needed
    plot.nominal.p="F",
    plot.FDR="F",
    
    random.seed=7392854
  ))



###Kaplan Meier curves and summary table
library(survival)

sig.col2=c(rgb(0,128,0,max=255),rgb(228,108,10,max=255))
sig.col3=c(rgb(0,128,0,max=255),rgb(127,127,127,max=255),rgb(228,108,10,max=255))

  cli <- read.delim(input.clinical.data, header=T)
  ntp<- read.delim(paste(ntp.output.name, "prediction_summary.txt",sep = "_"), header =T)
  d <- merge(cli, ntp, by.x =pt.ID, by.y ="sample.names")
  attach(d)
  d$risk.crude <- ifelse(d$predict.label ==1, "2_poor", "1_good")
  d$risk.p <- ifelse(d$predict.label.p0.05 ==1, "3_poor", ifelse(d$predict.label.p0.05 ==2, "1_good", "2_inte"))
  d$risk.p_poor <- ifelse(d$predict.label.p0.05 ==1, "2_poor", "1_rest")
  d$risk.p_good <- ifelse(d$predict.label.p0.05 ==2, "1_good", "2_rest")
  d$risk.FDR <- ifelse(d[, names(d)[grep("predict.label.FDR", names(d))]] ==1, "3_poor", ifelse(d[, names(d)[grep("predict.label.FDR", names(d))]] ==2, "1_good", "2_inte"))
  d$risk.FDR_poor <- ifelse(d[, names(d)[grep("predict.label.FDR", names(d))]] ==1, "2_poor", "1_rest")
  d$risk.FDR_good <- ifelse(d[, names(d)[grep("predict.label.FDR", names(d))]] ==2, "1_good", "2_rest")
  detach(d)
  eval(parse(text=paste("time.variable <- d$",time,sep="")))
  eval(parse(text=paste("censor.variable <- d$",censor,sep="")))
  
  
  ##KM curves 
  pdf(paste(paste(ntp.output.name, "KM.pdf",sep = "_")))
  par(mfrow=c(3,3))
  plot(survfit(Surv(time.variable, censor.variable)~risk.crude, data = d),ylim = ylim,conf.int = conf.int, mark.time = mark.time,fun=fun, col = sig.col2, main = "Crude")
  plot.new()
  plot.new()
  plot(survfit(Surv(time.variable, censor.variable)~risk.p, data = d),ylim = ylim,conf.int = conf.int, fun=fun, mark.time = mark.time, col = sig.col3, main = "p<0.05, 3 classes")
  plot(survfit(Surv(time.variable, censor.variable)~risk.p_poor, data = d),ylim = ylim,conf.int = conf.int, fun=fun, mark.time = mark.time, col = sig.col2, main = "p<0.05, poor vs the rest")
  plot(survfit(Surv(time.variable, censor.variable)~risk.p_good, data = d),ylim = ylim,conf.int = conf.int, fun=fun, mark.time = mark.time, col = sig.col2, main = "p<0.05, good vs the rest")
  plot(survfit(Surv(time.variable, censor.variable)~risk.FDR, data = d),ylim = ylim,conf.int = conf.int, fun=fun, mark.time = mark.time, col = sig.col3, main = "FDR, 3 classes")
  plot(survfit(Surv(time.variable, censor.variable)~risk.FDR_poor, data = d),ylim = ylim,conf.int = conf.int, fun=fun, mark.time = mark.time, col = sig.col2, main = "FDR, poor vs the rest")
  plot(survfit(Surv(time.variable, censor.variable)~risk.FDR_good, data = d),ylim = ylim,conf.int = conf.int, fun=fun, mark.time = mark.time, col = sig.col2, main = "FDR, good vs the rest")
  dev.off()
  
  #logrank_p

  if(length(unique(d$risk.crude))==1){crude.sdf <- NA}else{crude.sdf <- survdiff(Surv(time.variable, censor.variable)~risk.crude, data = d)}
  if(length(unique(d$risk.crude))==1){crude.p <- NA}  else{crude.p   <-  1 - pchisq(crude.sdf$chisq, length(crude.sdf$n) - 1)}
  
  if(length(unique(d$risk.p))==1){risk.p.sdf <- NA}else{risk.p.sdf <- survdiff(Surv(time.variable, censor.variable)~risk.p, data = d)}
  if(length(unique(d$risk.p))==1){risk.p.p <- NA}else{risk.p.p <-  1 - pchisq(risk.p.sdf$chisq, length(risk.p.sdf$n) - 1)}
  
  if(length(unique(d$risk.p_poor))==1){risk.p_poor.sdf <- NA}else{risk.p_poor.sdf <- survdiff(Surv(time.variable, censor.variable)~risk.p_poor, data = d)}
  if(length(unique(d$risk.p_poor))==1){risk.p_poor.p <- NA}else{risk.p_poor.p <-  1 - pchisq(risk.p_poor.sdf$chisq, length(risk.p_poor.sdf$n) - 1)}
  

  if(length(unique(d$risk.p_good))==1){risk.p_good.sdf <- NA}else{risk.p_good.sdf <- survdiff(Surv(time.variable, censor.variable)~risk.p_good, data = d)}
  if(length(unique(d$risk.p_good))==1){risk.p_good.p <- NA}else{risk.p_good.p <-  1 - pchisq(risk.p_good.sdf$chisq, length(risk.p_good.sdf$n) - 1)}
  
  
  if(length(unique(d$risk.FDR))==1){risk.FDR.sdf <- NA}else{risk.FDR.sdf <- survdiff(Surv(time.variable, censor.variable)~risk.FDR, data = d)}
  if(length(unique(d$risk.FDR))==1){risk.FDR.p <- NA}else{risk.FDR.p <-  1 - pchisq(risk.FDR.sdf$chisq, length(risk.FDR.sdf$n) - 1)}
  
  if(length(unique(d$risk.FDR_poor))==1){risk.FDR_poor.sdf <- NA}else{risk.FDR_poor.sdf <- survdiff(Surv(time.variable, censor.variable)~risk.FDR_poor, data = d)}
  if(length(unique(d$risk.FDR_poor))==1){risk.FDR_poor.p <- NA}else{risk.FDR_poor.p <-  1 - pchisq(risk.FDR_poor.sdf$chisq, length(risk.FDR_poor.sdf$n) - 1)}
  
  if(length(unique(d$risk.FDR_good))==1){risk.FDR_good.sdf <- NA}else{risk.FDR_good.sdf <- survdiff(Surv(time.variable, censor.variable)~risk.FDR_good, data = d)}
  if(length(unique(d$risk.FDR_good))==1){risk.FDR_good.p <- NA}else{risk.FDR_good.p <-  1 - pchisq(risk.FDR_good.sdf$chisq, length(risk.FDR_good.sdf$n) - 1)}

  logrank_pvalues <-  cox.summary <- matrix(c("NA","NA", "", "NA", "NA","NA", "", "NA", "NA"), nrow = 9, ncol =1)
    
  logrank_pvalues[1,] <- crude.p
  logrank_pvalues[2,] <- risk.p.p
                       
  logrank_pvalues[4,] <- risk.p_poor.p
  logrank_pvalues[5,] <- risk.p_good.p
  logrank_pvalues[6,] <- risk.FDR.p
                       
  logrank_pvalues[8,] <- risk.FDR_poor.p
  logrank_pvalues[9,] <- risk.FDR_good.p
                       
  
  ###Cox regression

  if(length(unique(d$risk.crude))==1){crude.cox <- rep(NA,5)}else{crude.cox  <- summary(coxph(Surv(time.variable, censor.variable)~risk.crude, data = d))$coefficients}
  if(length(unique(d$risk.p))==1){risk.p.cox <- matrix(rep(NA,10),nrow=2)}else{risk.p.cox  <- summary(coxph(Surv(time.variable, censor.variable)~risk.p, data = d))$coefficients}
  if(length(unique(d$risk.p_poor))==1){risk.p_poor.cox <- rep(NA,5)}else{risk.p_poor.cox  <- summary(coxph(Surv(time.variable, censor.variable)~risk.p_poor, data = d))$coefficients}
  if(length(unique(d$risk.p_good))==1){risk.p_good.cox <- rep(NA,5)}else{risk.p_good.cox  <- summary(coxph(Surv(time.variable, censor.variable)~risk.p_good, data = d))$coefficients}
  
  if(length(unique(d$risk.FDR))==1){risk.FDR.cox <- matrix(rep(NA,10),nrow=2)}else{risk.FDR.cox  <- summary(coxph(Surv(time.variable, censor.variable)~risk.FDR, data = d))$coefficients}
  if(length(unique(d$risk.FDR_poor))==1){risk.FDR_poor.cox <- rep(NA,5)}else{risk.FDR_poor.cox  <- summary(coxph(Surv(time.variable, censor.variable)~risk.FDR_poor, data = d))$coefficients}
  if(length(unique(d$risk.FDR_good))==1){risk.FDR_good.cox <- rep(NA,5)}else{risk.FDR_good.cox  <- summary(coxph(Surv(time.variable, censor.variable)~risk.FDR_good, data = d))$coefficients}
  
  
  crude.cox.summary          <- c(crude.cox[2],        exp(crude.cox[1]-1.96 * crude.cox[3]),          exp(crude.cox[1]+1.96 * crude.cox[3]), crude.cox[5])
  risk.p.cox_inter_summary   <- c(risk.p.cox[1,2],     exp(risk.p.cox[1,1]-1.96 * risk.p.cox[1,3]),    exp(risk.p.cox[1,1]+1.96 * risk.p.cox[1,3]), risk.p.cox[1, 5])
  risk.p.cox_poor_summary    <- c(risk.p.cox[2,2],     exp(risk.p.cox[2,1]-1.96 * risk.p.cox[2,3]),    exp(risk.p.cox[2,1]+1.96 * risk.p.cox[2,3]), risk.p.cox[2, 5])
  risk.p_poor.cox.summary    <- c(risk.p_poor.cox[2],  exp(risk.p_poor.cox[1]-1.96 * risk.p_poor.cox[3]), exp(risk.p_poor.cox[1]+1.96 * risk.p_poor.cox[3]), risk.p_poor.cox[5])
  risk.p_good.cox.summary    <- c(risk.p_good.cox[2],  exp(risk.p_good.cox[1]-1.96 * risk.p_good.cox[3]), exp(risk.p_good.cox[1]+1.96 * risk.p_good.cox[3]), risk.p_good.cox[5])
  risk.FDR.cox_inter_summary <- c(risk.FDR.cox[1,2],   exp(risk.FDR.cox[1,1]-1.96 * risk.FDR.cox[1,3]),   exp(risk.FDR.cox[1,1]+1.96 * risk.FDR.cox[1,3]), risk.FDR.cox[1, 5])
  risk.FDR.cox_poor_summary  <- c(risk.FDR.cox[2,2],    exp(risk.FDR.cox[2,1]-1.96 * risk.FDR.cox[2,3]),  exp(risk.FDR.cox[2,1]+1.96 * risk.FDR.cox[2,3]), risk.FDR.cox[2, 5])
  risk.FDR_poor.cox.summary  <- c(risk.FDR_poor.cox[2], exp(risk.FDR_poor.cox[1]-1.96 * risk.FDR_poor.cox[3]), exp(risk.FDR_poor.cox[1]+1.96 * risk.FDR_poor.cox[3]), risk.FDR_poor.cox[5])
  risk.FDR_good.cox.summary  <- c(risk.FDR_good.cox[2], exp(risk.FDR_good.cox[1]-1.96 * risk.FDR_good.cox[3]), exp(risk.FDR_good.cox[1]+1.96 * risk.FDR_good.cox[3]), risk.FDR_good.cox[5])
  
  cox.summary <- matrix(NA,nrow=9,ncol=4)
  
  cox.summary[1,] <- crude.cox.summary
  cox.summary[2,] <- risk.p.cox_inter_summary
  cox.summary[3,] <- risk.p.cox_poor_summary
  cox.summary[4,] <- risk.p_poor.cox.summary
  cox.summary[5,] <- risk.p_good.cox.summary
  cox.summary[6,] <- risk.FDR.cox_inter_summary
  cox.summary[7,] <- risk.FDR.cox_poor_summary
  cox.summary[8,] <- risk.FDR_poor.cox.summary
  cox.summary[9,] <- risk.FDR_good.cox.summary
                     
  row_names<- c("Crude",
                "p0.05_3cls_inter",
                "p0.05_3cls_poor", 
                "p0.05_3cls_poor_vs_rest",
                "p0.05_3cls_rest_vs_good", 
                "FDR_3cls_inter",
                "FDR_3cls_poor", 
                "FDR_3cls_poor_vs_rest",
                "FDR_3cls_rest_vs_good")
  
  num.cls <- c(as.numeric(table(d$risk.crude)[2]),
               as.numeric(table(d$risk.p)[2]), as.numeric(table(d$risk.p)[3]),
               as.numeric(table(d$risk.p_poor)[2]),
               as.numeric(table(d$risk.p_good)[2]),
               as.numeric(table(d$risk.FDR)[2]), as.numeric(table(d$risk.FDR)[3]),
               as.numeric(table(d$risk.FDR_poor)[2]),
               as.numeric(table(d$risk.FDR_good)[2]))
  
  d.summary <- cbind(row_names,num.cls, logrank_pvalues, cox.summary)
  colnames(d.summary) <- c("variables","N","Logrank_p", "HR", "HR_95%CI_low", "HR_95%CI_high", "Cox_p")
  
  
  
  write.table(d.summary, paste(ntp.output.name, "survival_analysis_summary.txt",sep = "_"), sep = "\t", row.names = F, col.names = T)
  save(d.summary, file=(paste(ntp.output.name, "survival_analysis_summary.RData",sep = "_")))
}


NTP<-function(
  # file I/O
  input.exp.filename,
  input.features.filename,
  output.name="NTP",
  
  norm.method="row.std", # "row.std.ref","ratio.ref"
  ref.sample.file=NULL,
  
  # temp.nn distance & row normalize
  distance.selection="cosine", # "correlation" or "cosine"
  weight.genes="T",   # only for 2 cls
  
  # resampling to generate null dist
  num.resamplings=1000,
  within.sig="F",
  
  # outputs
  GenePattern.output="T",
  signature.heatmap="T",
  p.sample.bar=0.05, # NA if not needed
  FDR.sample.bar=0.05, # NA if not needed
  plot.nominal.p="T",
  plot.FDR="T",
  
  random.seed=7392854
)
{
  
  #  suppressWarnings()
  
  # Advanced setting
  
  #  row.norm="T"
  col.range=3         # SD in heatmap
  heatmap.legend=signature.heatmap
  plot.distance="F"
  # histgram of null dist for the distance
  histgram.null.dist="F"
  hist.br=30
  
  # for dChip
  dchip.output="F"
  
  # for GenePattern
  num.resamplings<-as.numeric(num.resamplings)
  col.range<-as.numeric(col.range)
  hist.br<-as.numeric(hist.br)  # bin number for resampled dist histgram
  p.sample.bar <- as.numeric(p.sample.bar)
  FDR.sample.bar <- as.numeric(FDR.sample.bar)
  random.seed <- as.numeric(random.seed)
  #  if (FDR.sample.bar!="NA"){
  #    FDR.sample.bar <- as.numeric(FDR.sample.bar)
  #    if (is.numeric(FDR.sample.bar)==F){
  #      stop("### Provide numerical value (0~1) for FDR.sample.bar! ###")
  #    }
  #  }
  if (is.null(ref.sample.file)){
  }else{
    if (ref.sample.file=="NULL"){
      ref.sample.file <- NULL
    }
  }
  
  # set random seed
  set.seed(random.seed)
  
  ### input ###
  
  # selected features used for prediction
  features<-read.delim(input.features.filename,header=T,check.names=F)
  
  ## file format check
  if (length(features[1,])!=3 & length(features[1,])!=4){
    stop("### Please use features file format! ###")
  }
  if (length(features[1,])<4 & weight.genes=="T"){
    weight.genes <- "F"
  }
  third.col<-rownames(table(features[,3]))
  if (is.na(as.numeric(third.col[1]))){
    stop("### The 3rd column of feature file should be numerical! ###")
  }
  
  feat.col.names<-colnames(features)
  feat.col.names[1:2]<-c("ProbeID","GeneName")
  colnames(features)<-feat.col.names
  
  num.features<-length(features[,1])
  num.cls<-length(table(features[,3]))
  feature.col.num <- length(features[1,])
  
  ord<-seq(1:num.features)
  features<-cbind(ord,features)  # add order column to "features"
  
  # expression data
  ## file format check
  if (regexpr(".gct$",input.exp.filename)==-1){
    stop("### Gene expression data should be .gct format! ###")
  }
  exp.dataset<-read.delim(input.exp.filename,header=T,skip=2)
  #  exp.dataset<-read.delim(input.exp.filename,header=T,strip.white=T,dec=".",skip=2,check.names=F)
  colnames(exp.dataset)[1:2] <- c("ProbeID","GeneName")
  
  ## Other dataset's mean & SD for row normalization (optional)
  if (!is.null(ref.sample.file)){
    ref.sample <- read.delim(ref.sample.file,header=T)
    if (dim(ref.sample)[2]!=4 & is.numeric(ref.sample[1,3]) & is.numeric(ref.sample[1,4])){
      stop("### mean & SD file format incorrect! ###")
    }
    colnames(ref.sample)[1:4] <- c("ProbeID","SomeName","mean","sd")
    merged.dataset <- merge(ref.sample,exp.dataset,sort=F)
    
    ref.sample <- merged.dataset[,1:4]
    exp.dataset <- merged.dataset[,c(1,5:dim(merged.dataset)[2])]
  }
  
  ProbeID<-exp.dataset[,1]
  gene.names<-exp.dataset[,2]
  num.samples<-(length(exp.dataset[1,])-2)
  exp.dataset<-exp.dataset[-c(1:2)]
  
  exp.for.sample.names<-read.delim(input.exp.filename,header=F,skip=2)  # read sample names
  sample.names<-as.vector(as.matrix(exp.for.sample.names[1,3:length(exp.for.sample.names[1,])]))
  
  # row normalize
  
  normed.exp.dataset<-exp.dataset
  
  if (norm.method=="row.std"){
    exp.mean <- apply(exp.dataset,1,mean,na.rm=T)
    exp.sd <- apply(exp.dataset,1,sd,na.rm=T)
    normed.exp.dataset<-(exp.dataset-exp.mean)/exp.sd
  }
  if (norm.method=="row.std.ref"){
    if (is.null(ref.sample)){
      stop("### Provide reference sample data! ###")
    }
    exp.mean <- as.numeric(as.vector(ref.sample$mean))
    exp.sd <- as.numeric(as.vector(ref.sample$sd))
    normed.exp.dataset<-(exp.dataset-exp.mean)/exp.sd
  }
  if (norm.method=="ratio.ref"){
    if (is.null(ref.sample)){
      stop("### Provide reference sample data! ###")
    }
    exp.mean <- as.numeric(as.vector(ref.sample$mean))
    normed.exp.dataset<- exp.dataset/exp.mean
  }
  
  normed.exp.dataset<-cbind(ProbeID,normed.exp.dataset)
  
  # extract features from normed.exp.dataset
  
  exp.dataset.extract<-merge(features,normed.exp.dataset,sort=F)
  if (length(exp.dataset.extract[,1])<1){
    stop("### No matched probes! ###")
  }
  
  order.extract<-order(exp.dataset.extract[,2])
  exp.dataset.extract<-exp.dataset.extract[order.extract,]
  order.extract.after<-exp.dataset.extract[,2]
  exp.dataset.extract<-exp.dataset.extract[-2]
  
  if (weight.genes=="F"){
    features.extract<-exp.dataset.extract[,1:3]
    if (feature.col.num==4){
      exp.dataset.extract <- exp.dataset.extract[-4]
    }
    features.extract<-cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)
    num.features.extract<-length(features.extract[,1])
    
    ProbeID.extract<-as.vector(exp.dataset.extract[,1])
    exp.dataset.extract<-exp.dataset.extract[-c(1:3)]
    rownames(exp.dataset.extract)<-ProbeID.extract
  }
  
  #  weight.genes.vector <- rep(1,num.features)
  
  if (weight.genes=="T" & num.cls==2){
    features.extract<-exp.dataset.extract[,1:4]
    features.extract<-cbind(order.extract.after,features.extract) # order:ProbeID:gene name:cls:wt(if any)
    
    #    if (is.numeric(features[,4])){
    weight.genes.vector <- as.numeric(as.vector(features.extract[,5]))
    #    }else{
    if (is.numeric(weight.genes.vector)==F){
      stop("# Please use numeric values in 4th column!#")
    }
    
    num.features.extract<-length(features.extract[,1])
    
    ProbeID.extract<-as.vector(exp.dataset.extract[,1])
    exp.dataset.extract<-exp.dataset.extract[-c(1:4)]
    rownames(exp.dataset.extract)<-ProbeID.extract
  }
  
  # make template
  
  for (i in 1:num.cls){
    temp.temp<-as.numeric(as.vector(features.extract[,4]))
    temp.temp[temp.temp!=i]<-0
    temp.temp[temp.temp==i]<-1
    eval(parse(text=paste("temp.",i,"<-temp.temp",sep="")))
    #    eval(parse(text=paste("temp\.",i,"<-temp\.temp",sep="")))  ### for < R-2.4.0
  }
  
  # weighted template (only for 2cls)
  
  if (weight.genes=="T" & num.cls==2){
    temp.1 <- weight.genes.vector
    temp.2 <- -weight.genes.vector
  }
  
  ### compute distance and p-value ###
  
  predict.label<-vector(length=num.samples,mode="numeric")
  dist.to.template<-vector(length=num.samples,mode="numeric")
  dist.to.cls1<-vector(length=num.samples,mode="numeric")
  
  rnd.feature.matrix<-matrix(0,nrow=num.features.extract,ncol=num.resamplings)
  
  perm.dist.vector<-vector(length=num.resamplings*num.cls,mode="numeric")
  nominal.p<-vector(length=num.samples,mode="numeric")
  BH.FDR<-vector(length=num.samples,mode="numeric")
  Bonferroni.p<-vector(length=num.samples,mode="numeric")
  
  for (i in 1:num.samples){
    
    print(paste("sample # ",i,sep=""))
    
    current.sample <- as.vector(exp.dataset.extract[,i])
    
    # compute original distance
    
    orig.dist.to.all.temp <- vector(length=num.cls,mode="numeric")
    
    if (weight.genes=="T"){   # weight sample data
      current.sample <- current.sample*abs(weight.genes.vector)
    }
    
    if (distance.selection=="cosine"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        #        eval(parse(text=paste("current\.temp <- temp\.",o,sep="")))  ### for < R-2.4.0
        orig.dist.to.all.temp[o]<-sum(current.temp*current.sample)/
          (sqrt(sum(current.temp^2))*sqrt(sum(current.sample^2)))
      }
    }
    if (distance.selection=="correlation"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        #        eval(parse(text=paste("current\.temp <- temp\.",o,sep="")))  ### for < R-2.4.0
        orig.dist.to.all.temp[o] <- cor(current.temp,current.sample,method="pearson",use="complete.obs")
      }
    }
    
    if (num.cls==2){           # find nearest neighbor (2 classes)
      if (orig.dist.to.all.temp[1]>=orig.dist.to.all.temp[2]){
        predict.label[i]<-1
        dist.to.template[i]<-1-orig.dist.to.all.temp[1]
        dist.to.cls1[i]<--(orig.dist.to.all.temp[1]+1)
      }
      if (orig.dist.to.all.temp[1]<orig.dist.to.all.temp[2]){
        predict.label[i]<-2
        dist.to.template[i]<-1-orig.dist.to.all.temp[2]
        dist.to.cls1[i]<-orig.dist.to.all.temp[2]+1
      }
    }
    
    if (num.cls>2){
      for (o in 1:num.cls){       # find nearest neighbor (>2 classes)
        if (is.na(orig.dist.to.all.temp[o])!=T){
          if (orig.dist.to.all.temp[o]==max(orig.dist.to.all.temp,na.rm=T)){
            predict.label[i]<-o
            dist.to.template[i]<-1-orig.dist.to.all.temp[o]
            dist.to.cls1[i]<-(1-orig.dist.to.all.temp[o])+o
          }
        }
      }
    }
    
    # permutation test
    
    if (within.sig=="F"){     # generate resampled features from all probes
      for (p in 1:num.resamplings){
        rnd.feature.matrix[,p]<-sample(normed.exp.dataset[,(i+1)],num.features.extract,replace=F)
      }
    }
    if (within.sig=="T"){     # generate resampled features from only signature genes
      for (p in 1:num.resamplings){
        rnd.feature.matrix[,p]<-sample(exp.dataset.extract[,i],num.features.extract,replace=F)
      }
    }
    
    if (weight.genes=="T" & num.cls==2){
      rnd.feature.matrix <- rnd.feature.matrix*abs(weight.genes.vector)
    }
    
    # compute distance to all templates
    if (distance.selection=="cosine"){          # cosine
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        
        prod.sum<-apply(t(t(rnd.feature.matrix)*temp.resmpl),2,sum)
        
        data.sq.sum<-apply(rnd.feature.matrix^2,2,sum)
        temp.sq.sum<-sum(temp.resmpl^2)
        
        perm.dist.vector[(1+(num.resamplings*(res-1))):(num.resamplings*res)]<-
          (1-(prod.sum/(sqrt(data.sq.sum)*sqrt(temp.sq.sum))))
      }
    }
    
    if (distance.selection=="correlation"){          # correlation
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        perm.dist.vector[(1+(num.resamplings*(res-1))):(num.resamplings*res)]<-
          (1-as.vector(cor(rnd.feature.matrix,temp.resmpl,method="pearson",use="complete.obs")))
      }
    }
    
    # compute nominal p-value
    
    combined.stats.rank<-rank(c(dist.to.template[i],perm.dist.vector))
    nominal.p[i]<-combined.stats.rank[1]/length(combined.stats.rank)
    
    # histgram of combined null distributions
    
    if (histgram.null.dist=="T" & capabilities("png")==T){
      png(paste("resampled_",distance.selection,"_dist_histgram_",sample.names[i],".png",sep=""), type="cairo")
      hist(c(dist.to.template[i],perm.dist.vector),br=hist.br,main=paste(sample.names[i],", # resampling: ",num.resamplings,sep=""))
      dev.off()
    }
    
  } # main sample loop END
  
  # MCT correction
  
  BH.FDR<-nominal.p*num.samples/rank(nominal.p)
  Bonferroni.p<-nominal.p*num.samples
  
  BH.FDR[BH.FDR>1]<-1
  Bonferroni.p[Bonferroni.p>1]<-1
  
  # prediction results w/ prediction confidence
  
  predict.label.p <- predict.label.FDR <- predict.label
  predict.label.p[which(nominal.p>=p.sample.bar)] <- 0
  predict.label.FDR[which(BH.FDR>=FDR.sample.bar)] <- 0
  
  ### output ###
  
  # prediction results
  
  dist.to.cls1.rank <- rank(dist.to.cls1)
  pred.summary <- cbind(sample.names,predict.label,predict.label.p,predict.label.FDR,dist.to.template,dist.to.cls1.rank,
                        nominal.p,BH.FDR,Bonferroni.p)
  colnames(pred.summary)[3] <- paste("predict.label.p",p.sample.bar,sep="")
  colnames(pred.summary)[4] <- paste("predict.label.FDR",FDR.sample.bar,sep="")
  
  write.table(pred.summary,paste(output.name,"_prediction_summary.txt",sep=""),
              quote=F,sep="\t",row.names=F)
  
  # extracted features
  
  if (weight.genes=="T" & num.cls==2){
    write.table(features.extract[,2:5],paste(output.name,"_features.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  if (weight.genes=="F"){
    write.table(features.extract[,2:4],paste(output.name,"_features.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  
  # sorted exp dataset for heatmap (row normalized)
  
  t.dataset<-t(exp.dataset.extract)               # sort samples
  t.dataset<-cbind(dist.to.cls1,t.dataset)
  ts.dataset<-t.dataset[order(t.dataset[,1]),]
  to.dataset.out<-ts.dataset[,2:(num.features.extract+1)]
  
  sorted.dataset<-t(to.dataset.out)
  heatmap.dataset<-as.matrix(sorted.dataset)
  #  if (.Platform$OS.type == "windows") {
  #    sorted.dataset<-matrix(gsub(" ","",sorted.dataset),ncol=(num.samples+2))
  #  }
  
  # sorted exp dataset for spreadsheets (not normalized)
  
  exp.dataset.gannot<-cbind(ProbeID,exp.dataset)
  exp.dataset.extract<-merge(features.extract,exp.dataset.gannot,sort=F) # redefine exp.dataset.extract
  
  order.extract<-order(exp.dataset.extract[,2])   # sort genes
  exp.dataset.extract<-exp.dataset.extract[order.extract,]
  if (weight.genes=="T" & num.cls==2){
    exp.dataset.extract<-exp.dataset.extract[-c(1:5)]
  }
  if (weight.genes=="F"){
    exp.dataset.extract<-exp.dataset.extract[-c(1:4)]
  }
  
  t.dataset<-t(exp.dataset.extract)               # sort samples
  t.dataset<-cbind(dist.to.cls1,t.dataset)
  ts.dataset<-t.dataset[order(t.dataset[,1]),]
  to.dataset.out<-ts.dataset[,2:(num.features.extract+1)]
  
  sorted.dataset<-t(to.dataset.out)
  sorted.dataset<-cbind(features.extract[,2:3],sorted.dataset)
  
  sorted.dataset.header<-c("ProbeID","GeneName",sample.names[order(t.dataset[,1])])
  sorted.dataset<-t(cbind(sorted.dataset.header,t(sorted.dataset)))
  
  if (.Platform$OS.type == "windows") {
    sorted.dataset<-matrix(gsub(" ","",sorted.dataset),ncol=(num.samples+2))
  }
  
  # output for GenePattern
  
  if (GenePattern.output=="T"){
    
    # exp data
    write.table("#1.2",paste(output.name,"_sorted.dataset.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F)
    write.table(paste(num.features.extract,num.samples,sep="\t"),paste(output.name,"_sorted.dataset.gct",sep="")
                ,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table(sorted.dataset,paste(output.name,"_sorted.dataset.gct",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F,append=T)
    
    # cls ffiles (unsorted, sorted)
    cls.out<-matrix(0,nrow=3,ncol=1)
    
    ## unsorted cls
    cls.out[1,]<-paste(num.samples," ",num.cls," 1",sep="")  # line 1
    cls.out[2,]<-paste("# ",paste(unique(predict.label),collapse=" ")) # line 2
    predict.label.out<-as.numeric(predict.label)-1                   # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")
    
    write.table(cls.out,paste(output.name,"_predicted_unsorted.cls",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    ## unsorted cls (with non-confident prediction by nominal p)
    cls.out[1,]<-paste(num.samples," ",(num.cls+1)," 1",sep="")  # line 1
    cls.out[2,]<-paste("# ",paste(unique(predict.label.p),collapse=" ")) # line 2
    predict.label.out<-as.numeric(predict.label.p)                     # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")
    
    write.table(cls.out,paste(output.name,"_predicted_unsorted_p",p.sample.bar,".cls",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    ## unsorted cls (with non-confident prediction by FDR)
    cls.out[1,]<-paste(num.samples," ",(num.cls+1)," 1",sep="")  # line 1
    cls.out[2,]<-paste("# ",paste(unique(predict.label.FDR),collapse=" ")) # line 2
    predict.label.out<-as.numeric(predict.label.FDR)                     # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")
    
    write.table(cls.out,paste(output.name,"_predicted_unsorted_FDR",FDR.sample.bar,".cls",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    ## sorted cls
    sorted.predict.label<-sort(as.numeric(predict.label))
    cls.out[2,]<-paste("# ",paste(unique(sorted.predict.label),collapse=" ")) # line 2
    predict.label.out<-sorted.predict.label-1    # line 3
    cls.out[3,]<-paste(predict.label.out,collapse=" ")
    
    write.table(cls.out,paste(output.name,"_predicted_sorted.cls",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
  }
  
  # output for dChip
  
  if (dchip.output=="T"){
    
    # exp data
    write.table(sorted.dataset,paste(output.name,"_dChip_sorted.dataset.txt",sep=""),
                quote=F,sep="\t",row.names=F,col.names=F)
    
    # sample info
    sample.info<-cbind(sample.names,sample.names,predict.label)
    write.table(sample.info,paste(output.name,"_dChip_sample_info.txt",sep=""),
                quote=F,sep="\t",row.names=F)
    
    # gene info
    dchip.gene.info<-cbind(features.extract[,2:3],NA,features.extract[,3],features.extract[,4],NA)
    colnames(dchip.gene.info)<-c("Probe Set","Identifier","FirstOfLocuslink","FirstOfName","FirstOfFunction","Description")
    write.table(dchip.gene.info,paste(output.name,"_dChip_gene_info.txt",sep=""),
                quote=F,sep="\t",row.names=F)
  }
  
  # heatmap
  
  if (signature.heatmap=="T" & capabilities("png")==T){
    
    subclass.col.source <- c("red","blue","yellow","green","purple","orange","lightblue","darkgreen")
    predict.col.vector <- unique(sort(predict.label))
    subclass.col <- subclass.col.source[predict.col.vector]
    
    heatmap.col <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000")
    
    #    if (row.norm=="T"){
    #      exp.mean <- apply(heatmap.dataset,1,mean) #,na.rm=T)
    #      exp.sd <- apply(heatmap.dataset,1,sd) #,na.rm=T)
    #      heatmap.dataset <- (heatmap.dataset-exp.mean)/exp.sd
    heatmap.dataset[heatmap.dataset>col.range] <- col.range
    heatmap.dataset[heatmap.dataset< -col.range] <- -col.range  # bug fixed 06/18/2012
    #    }
    
    ncol.heat <- length(heatmap.dataset[1,])
    nrow.heat <- length(heatmap.dataset[,1])
    
    heatmap.dataset <- apply(heatmap.dataset,2,rev)
    
    num.pred <- as.vector(table(predict.label))
    num.pred.gene <- as.vector(table(features.extract[,4]))
    
    increment.sample <- cumsum(num.pred)
    increment.sample <- c(0,increment.sample)
    
    increment.gene <- cumsum(num.pred.gene)
    increment.gene <- c(0,increment.gene)
    
    png(paste(output.name,"_heatmap.png",sep=""), type="cairo")
    image(1:ncol.heat,1:nrow.heat,t(heatmap.dataset),axes=F,col=heatmap.col,zlim=c(-col.range,col.range),xlim=c(-0.5,(ncol.heat+0.5+round(ncol.heat*0.05))),ylim=c(-0.5,(nrow.heat+0.5+round(nrow.heat*0.08))),xlab=NA,ylab=NA)
    
    for (c in 1:num.cls){                          # gene bar
      rect((ncol.heat+1),0.5,(ncol.heat+0.5+round(ncol.heat*0.05)),(nrow.heat+0.5-increment.gene[c]),col=subclass.col[c],xpd=T,border=F)
    }
    for (c in 1:length(num.pred)){               # sample bar
      rect((0.5+increment.sample[c]),(nrow.heat+2),(ncol.heat+0.5),(nrow.heat+0.5+round(nrow.heat*0.08)),col=subclass.col[c],xpd=T,border=F)
    }
    dev.off()
    
    # heatmap legend
    
    if (heatmap.legend=="T"){
      png(paste(output.name,"_heatmap_legend.png",sep=""), type="cairo")
      par(plt=c(.1,.9,.45,.5))
      a=matrix(seq(1:15),nc=1)
      image(a,col=heatmap.col,xlim=c(0,1),axes=F,yaxt="n")
      box()
      dev.off()
    }
    
    # p-value sample bar
    
    #    if (p.sample.bar!="NA" & num.cls==2){
    if (num.cls==2){
      p.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      p.bar.vector[which(nominal.p[order(dist.to.cls1.rank)]>=p.sample.bar)] <- 3
      png(paste(output.name,"_p_",p.sample.bar,"_sample_bar.png",sep=""), type="cairo")
      par(plt=c(.1,.9,.45,.5))
      a=matrix(p.bar.vector,nc=1)
      #        image(a,col=c("red","blue","gray"),xlim=c(0,1),axes=F,yaxt="n")
      image(a,col=c("red","blue","gray"),axes=F,yaxt="n")
      dev.off()
    }
    
    #    if (p.sample.bar!="NA" & num.cls>2){
    if (num.cls>2){
      p.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      p.bar.vector[which(nominal.p[order(dist.to.cls1.rank)]>=p.sample.bar)] <- (num.cls+1)
      uni.p.bar.vector <- sort(unique(p.bar.vector))
      if (length(uni.p.bar.vector)>1){
        n.sig.cls <- length(uni.p.bar.vector)-1
        uni.p.bar.vector <- uni.p.bar.vector[1:n.sig.cls]
      }else{
        uni.p.bar.vector <- NULL
      }
      
      sig.subclass.col <- c(subclass.col.source[1:num.cls],"gray")
      png(paste(output.name,"_p_",p.sample.bar,"_sample_bar.png",sep=""), type="cairo")
      par(plt=c(.1,.9,.45,.5))
      a=matrix(p.bar.vector,nc=1)
      image(a,col=sig.subclass.col,axes=F,yaxt="n")
      dev.off()
    }
    
    # FDR sample bar
    
    #    if (FDR.sample.bar!="NA" & num.cls==2){
    if (num.cls==2){
      fdr.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      fdr.bar.vector[which(BH.FDR[order(dist.to.cls1.rank)]>=FDR.sample.bar)] <- 3
      png(paste(output.name,"_FDR_",FDR.sample.bar,"_sample_bar.png",sep=""), type="cairo")
      par(plt=c(.1,.9,.45,.5))
      a=matrix(fdr.bar.vector,nc=1)
      #        image(a,col=c("red","blue","gray"),xlim=c(0,1),axes=F,yaxt="n")
      image(a,col=c("red","blue","gray"),axes=F,yaxt="n")
      dev.off()
    }
    
    #    if (FDR.sample.bar!="NA" & num.cls>2){
    if (num.cls>2){
      fdr.bar.vector <- predict.label[order(dist.to.cls1.rank)]
      fdr.bar.vector[which(BH.FDR[order(dist.to.cls1.rank)]>=FDR.sample.bar)] <- (num.cls+1)
      uni.fdr.bar.vector <- sort(unique(fdr.bar.vector))
      if (length(uni.fdr.bar.vector)>1){
        n.sig.cls <- length(uni.fdr.bar.vector)-1
        uni.fdr.bar.vector <- uni.fdr.bar.vector[1:n.sig.cls]
      }else{
        uni.fdr.bar.vector <- NULL
      }
      
      sig.subclass.col <- c(subclass.col.source[1:num.cls],"gray")
      png(paste(output.name,"_FDR_",FDR.sample.bar,"_sample_bar.png",sep=""), type="cairo")
      par(plt=c(.1,.9,.45,.5))
      a=matrix(fdr.bar.vector,nc=1)
      image(a,col=sig.subclass.col,axes=F,yaxt="n")
      dev.off()
    }
    
  }
  
  # plot FDR
  
  if (plot.FDR=="T" & capabilities("png")==T){
    png(paste(output.name,"_FDR.png",sep=""), type="cairo")
    par(plt=c(0.1,0.95,0.4,0.6),las=2)
    plot(BH.FDR[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main="BH-FDR")
    box(lwd=2)
    dev.off()
  }
  
  # plot nominal-p
  
  if (plot.nominal.p=="T" & capabilities("png")==T){
    png(paste(output.name,"_nominal.p.png",sep=""), type="cairo")
    par(plt=c(0.1,0.95,0.4,0.6),las=2)
    plot(nominal.p[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main="nominal p-value")
    box(lwd=2)
    dev.off()
  }
  
  # plot distance to template
  
  if (plot.distance=="T" & capabilities("png")==T){
    png(paste(output.name,"_distance.png",sep=""), type="cairo")
    par(plt=c(0.1,0.95,0.4,0.6),las=2)
    plot(dist.to.template[order(dist.to.cls1)],pch=3,col="blue",lwd=2,ylim=c(0,1),main=paste("1 - ",distance.selection,sep=""))
    box(lwd=2)
    dev.off()
  }
  
}  # END main
