##############################################################
# SurvivalGene.R (ver.2)                  Aug.13,2008
#     
#     input:  (1) gct
#             (2) clinical data
#                *.txt, 1st column is sample name in the same
#                 order with .gct
#                *
#     output: (1) result
#             (2) mean & sd of empirical dist for the statistic
#                 (used for loocv.nn.R & nn.surv.iter.Kfold.R)
#                 Shapiro-Wilk normality test (W & p-value)
##############################################################

SurvivalGene <- function(
  input.filename.gct,
  input.filename.clinical,
    check.names.gct="F",
  output.file="SurvivalGene",

  time.field="time",
  censor.field="status",

  statistic.selection="cox.score", #  "cox.score", "coxph"
  trim.percent.2.side=0,
  nperm=1000,
  rnd.seed=56438219,

  emp.stat.dist="T"  # output mean & sd of empirical distribution of statistic
  )
{

  # for GenePattern

  trim.percent.2.side<-as.numeric(trim.percent.2.side)
  rnd.seed <- as.numeric(rnd.seed)

  set.seed(rnd.seed)

  # read input data

  if (check.names.gct=="F"){
    exp.all <- read.delim(input.filename.gct,header=T,skip=2,check.names=F)
  }else{
    exp.all <- read.delim(input.filename.gct,header=T,skip=2)
  }

  ## file format check
  if (regexpr(".gct$",input.filename.gct)==-1){
    stop("### Gene expression data should be .gct format! ###")
  }

  probeid <- as.vector(exp.all[,1])
  gene.names <- as.vector(exp.all[,2])
  num.samples<-(length(exp.all[1,])-2)
  num.genes <- length(exp.all[,1])
  exp.all<-exp.all[-c(1:2)]
  sample.names.all <- colnames(exp.all)
  rownames(exp.all) <- probeid

  clin.all <- read.delim(input.filename.clinical,header=T)

  eval(parse(text=paste("time <- clin.all$",time.field,sep="")))
  eval(parse(text=paste("censor <- clin.all$",censor.field,sep="")))

  sample.names.clin <- as.vector(clin.all[,1])

  ## check sample order btwn exp & clin
  for (i in 1:num.samples){
    if (sample.names.all[i]!=sample.names.clin[i]){
      stop("### Sample names don't match! ###")
    }
  }

  # trim survival time outliers

  num.trim.1.side <- round((num.samples*trim.percent.2.side+0.5)/2)
  order.survival <- order(time)
  trimmed.ordered.index <- order.survival[(num.trim.1.side+1):(num.samples-num.trim.1.side)]
  trimmed.time <- time[trimmed.ordered.index]
  trimmed.censor <- censor[trimmed.ordered.index]
  trimmed.exp <- exp.all[,trimmed.ordered.index]
  num.samples.trimmed <- length(trimmed.time)

  # set permutation matrices

  perm.stat <- matrix(0,nrow=num.genes,ncol=(nperm+1))
  perm.index <- matrix(0,nrow=num.samples.trimmed,ncol=nperm)
  perm.index.vector <- seq(1:num.samples.trimmed)
  for (i in 1:nperm){
    perm.index[,i] <- sample(perm.index.vector)
  }

  # cases at risk & death
  event.index <- which(trimmed.censor==1)
  event.time <- unique(trimmed.time[event.index])
  num.event.time <- length(event.time)

  # compute statistic
  if (statistic.selection=="cox.score"){

    # original statistics
    numerator <- 0
    denominator <- 0

    for (t in 1:num.event.time){
      index.at.risk <- which(trimmed.time>=event.time[t])
      num.at.risk <- length(index.at.risk)
      index.death <- which(trimmed.time==event.time[t])
      num.death <- length(index.death)

      if (is.null(dim(trimmed.exp[,index.at.risk]))==T){
        if (num.death==1){
          numerator <- numerator + trimmed.exp[,index.death]-num.death*trimmed.exp[,index.at.risk]
        }
       # if (num.death>1){
       #   numerator <- numerator + apply(trimmed.exp[,index.death],1,sum)-num.death*apply(trimmed.exp[,index.at.risk],1,mean)
       # }
        denominator <- denominator + (num.death/num.at.risk)*as.vector((trimmed.exp[,index.at.risk]-as.vector(trimmed.exp[,index.at.risk]))^2)
      }else{
        if (num.death==1){
          numerator <- numerator + trimmed.exp[,index.death]-num.death*apply(trimmed.exp[,index.at.risk],1,mean)
        }
        if (num.death>1){
          numerator <- numerator + apply(trimmed.exp[,index.death],1,sum)-num.death*apply(trimmed.exp[,index.at.risk],1,mean)
        }
        denominator <- denominator + (num.death/num.at.risk)*apply((trimmed.exp[,index.at.risk]-apply(trimmed.exp[,index.at.risk],1,mean))^2,1,sum)
#        denominator <- denominator + (num.death/num.at.risk)*as.vector(apply((trimmed.exp[,index.at.risk]-as.vector(apply(trimmed.exp[,index.at.risk],1,mean)))^2,1,sum))
      }      
    }
    perm.stat[,1] <-numerator/sqrt(denominator)

    # permutation startistic

    for (p in 1:nperm){

      print(paste("# Start permutation ",p,sep=""))

      perm.trimmed.exp <- trimmed.exp[,perm.index[,p]]

      numerator <- 0
      denominator <- 0

      for (t in 1:num.event.time){
        index.at.risk <- which(trimmed.time>=event.time[t])
        num.at.risk <- length(index.at.risk)
        index.death <- which(trimmed.time==event.time[t])
        num.death <- length(index.death)

        if (is.null(dim(perm.trimmed.exp[,index.at.risk]))==T){
          if (num.death==1){
            numerator <- numerator + perm.trimmed.exp[,index.death]-num.death*perm.trimmed.exp[,index.at.risk]
          }
         # if (num.death>1){
         #   numerator <- numerator + apply(perm.trimmed.exp[,index.death],1,sum)-num.death*apply(perm.trimmed.exp[,index.at.risk],1,mean)
         # }
          denominator <- denominator + (num.death/num.at.risk)*(perm.trimmed.exp[,index.at.risk]-perm.trimmed.exp[,index.at.risk])^2
        }else{
          if (num.death==1){
            numerator <- numerator + perm.trimmed.exp[,index.death]-num.death*apply(perm.trimmed.exp[,index.at.risk],1,mean)
          }
          if (num.death>1){
            numerator <- numerator + apply(perm.trimmed.exp[,index.death],1,sum)-num.death*apply(perm.trimmed.exp[,index.at.risk],1,mean)
          }
          denominator <- denominator + (num.death/num.at.risk)*apply((perm.trimmed.exp[,index.at.risk]-apply(perm.trimmed.exp[,index.at.risk],1,mean))^2,1,sum)
        }
      }
      perm.stat[,(p+1)] <-numerator/sqrt(denominator)
    }   # permutation statistics END

    # p-value

    rank.perm.stat <- nperm+1-t(apply(abs(perm.stat),1,rank,ties.method="min"))
    nominal.p <- rank.perm.stat[,1]/nperm

    rank.nom.p <- rank(nominal.p,ties.method="min")
    BH.FDR <- nominal.p*num.genes/rank.nom.p
    BH.FDR[BH.FDR>1] <- 1
    FWER <- nominal.p*num.genes
    FWER[FWER>1] <- 1

  } # cox.score END

  if (statistic.selection=="coxph"){
    library(survival)
    statistic <- p.value <- vector(length=num.genes,mode="numeric")
    for (g in 1:num.genes){
      print(paste("# Start gene ",g,"/",num.genes,sep=""))
      coxfit <- coxph(Surv(trimmed.time,trimmed.censor)~as.vector(as.matrix(trimmed.exp[g,])))
      statistic[g] <- as.numeric(coxfit$coef)
    }

  }  # coxph END

  # output
  if (statistic.selection=="cox.score"){
    statistic <- perm.stat[,1]
    output <- cbind(probeid,gene.names,statistic,nominal.p,BH.FDR,FWER)
  }
  if (statistic.selection=="coxph"){
    output <- cbind(probeid,gene.names,statistic,p.value)
  }
  write.table(output,paste(output.file,".txt",sep=""),quote=F,sep="\t",row.names=F)

  if (emp.stat.dist=="T"){
    mean <- apply(perm.stat[,2:(nperm+1)],1,mean)
    sd <- apply(perm.stat[,2:(nperm+1)],1,sd)
#    mean <- apply(perm.stat[,2:(nperm+1)],1,mean,na.rm=T)
#    sd <- apply(perm.stat[,2:(nperm+1)],1,sd,na.rm=T)
    sk.test.W <- sk.test.p <- vector(length=num.genes,mode="numeric")
    for (i in 1:num.genes){
      sk.test <- shapiro.test(as.vector(perm.stat[i,2:(nperm+1)]))
      sk.test.W[i] <- sk.test$statistic
      sk.test.p[i] <- sk.test$p.value
    }

    emp.stat.output <- cbind(probeid,mean,sd,sk.test.W,sk.test.p)
    write.table(emp.stat.output,paste(output.file,"_emp.stat.txt",sep=""),quote=F,sep="\t",row.names=F)
  }

#  return(list(output,perm.stat))

} # main END




  



