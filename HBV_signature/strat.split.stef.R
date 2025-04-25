##############################################################
# strat.split.stef.R                        4/1/2007
#   Added: L95,L96 for number of samples in split sets, 5/1/2014
#   Corretion: L77:"!="->"==",L172:"strat.param)", 5/1/2014
#              L100:"num.samples"->"num.exp.train.samples" 5/1/2014
#              L103:"num.sample"->"num.exp.test.samples" 5/1/2014
#     input:  gct
#             clinical data (.txt, 1st column is sample name
# currently only 2 splits
##############################################################

strat.split.stef <- function(
  # file I/O
  input.filename.gct,
  input.filename.clinical,

  output.name="Split",

  # split Kgrp parameters
  kfold=2,
  cls.var = "k.death",
  strata.var="SURVIVAL",
  strata.param=10,
  method=1,

  num.split=5,

  rnd.seed=2006
  )
{

  # for GenePattern

  kfold <- as.numeric(kfold)
  strata.param <- as.numeric(strata.param)
  method <- as.numeric(method)
  rnd.seed <- as.numeric(rnd.seed)

  # read input data

  exp.all <- read.delim(input.filename.gct,header=T,skip=2,check.names=F)

  ## file format check
  if (regexpr(".gct$",input.filename.gct)==-1){
    stop("### Gene expression data should be .gct format! ###")
  }

  num.genes <- length(exp.all[,1])
  num.samples <- length(exp.all[1,])-2

  sample.names.all <- colnames(exp.all)[3:length(exp.all[1,])]

  clin.all <- read.delim(input.filename.clinical,header=T)
  sample.names.clin <- as.vector(clin.all[,1])

  ## check sample order btwn exp & clin
  for (i in 1:num.samples){
    if (sample.names.all[i]!=sample.names.clin[i]){
      stop("### Sample names don't match! ###")
    }
  }

  set.seed(rnd.seed)

  for (spl in 1:num.split){
    print(paste("# Start split ",spl,sep=""))

    # split into k-folds

    clin.all.Kgrp <- split.Kgrp(clin.all,kfold=kfold,extreme.var=cls.var,
                            strata.var=strata.var,strata.param=strata.param, method=method)

    # assign samples to training or test set

    Kgrp.assigned <- clin.all.Kgrp$Kgrp

    index.train <- which(Kgrp.assigned==1)
    index.test <- which(Kgrp.assigned==2)

    exp.train <- exp.all[,c(1,2,(index.train+2))]
    exp.test <- exp.all[,c(1,2,(index.test+2))]

    clin.train <- clin.all[index.train,]
    clin.test <- clin.all[index.test,]

#    eval(parse(text=paste("sample.cls.train <- as.vector(clin.train$",cls.var,")",sep="")))
#    eval(parse(text=paste("sample.cls.test <- as.vector(clin.test$",cls.var,")",sep="")))

    exp.train.name <- paste(output.name,"_",spl,"_train.gct",sep="")
    exp.test.name <- paste(output.name,"_",spl,"_test.gct",sep="")

    clin.train.name <- paste(output.name,"_",spl,"_train_clinical.txt",sep="")
    clin.test.name <- paste(output.name,"_",spl,"_test_clinical.txt",sep="")

    num.exp.train.samples <- length(exp.train[1,])-2
    num.exp.test.samples <- length(exp.test[1,])-2

    old.warn<-options("warn"=-1)
    write.table("#1.2",exp.train.name,quote=F,sep="\t",row.names=F,col.names=F)
    write.table(paste(num.genes,num.exp.train.samples,sep="\t"),exp.train.name,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table(exp.train,exp.train.name,quote=F,sep="\t",row.names=F,col.names=T,append=T)

    write.table("#1.2",exp.test.name,quote=F,sep="\t",row.names=F,col.names=F)
    write.table(paste(num.genes,num.exp.test.samples,sep="\t"),exp.test.name,quote=F,sep="\t",row.names=F,col.names=F,append=T)
    write.table(exp.test,exp.test.name,quote=F,sep="\t",row.names=F,col.names=T,append=T)
    options(old.warn)

    write.table(clin.train,clin.train.name,quote=F,sep="\t",row.names=F,col.names=T)
    write.table(clin.test,clin.test.name,quote=F,sep="\t",row.names=F,col.names=T)
  }
}

split.Kgrp <- function(object, kfold=10, extreme.var = "extreme", strata.var="fup.month"
                       ,strata.param=10, method=1)
  {
    object$ind = 1:nrow(object)
    strataVar = NULL

    if(!is.null(strata.var))
      strataVar = object[,strata.var]

    ## split by extreme status: sampling is performed within Indolent and Lethal

    obj.split <- split(object,object[,extreme.var])

    ## method = 1

    setKgrp <- function(x)
      {
        ll = length(x$ind)

        if(!ll) ##length == 0
          return(NULL)

        rll = ll %% kfold

        ind = sample(x$ind,size=length(x$ind)) ## scramble ids
        ## number of samples multiple of kgrp
        if(!rll)
          kgrp = cbind(ind, rep(1:kfold,ll %/% kfold))
        ## number of samples not multiple of kgrp
        else
          {

            rgrp = sample(1:kfold,size=rll,replace=TRUE)
            agrp = c(rep(1:kfold,ll %/% kfold),rgrp)
            kgrp = cbind(ind,agrp)
          }

        return(kgrp)
      }

   ## Method = 2

    setKgrp2 <- function(x)
      {
        ll = length(x$ind)
        ## every subject randomly gets a number between 1:kfold, i.e
        ## which grp it'll belong to
        kgrp = cbind(x$ind, sample(1:kfold,ll,replace=TRUE))
        return(kgrp)
      }


    if(!is.null(strata.var))
      {
        if(is.null(strata.param))
          strata.param = 10

        ## compute breaks
        if(length(strata.param) == 1) ## number of quantiles
          split.breaks <- quantile(strataVar,seq(0,1,1/strata.param), type = 8)
        else
          split.breaks = strata.param



        sampleStrata = function(x)
          {
            strato = cut(x[,strata.var],breaks=split.breaks,right=TRUE,include=TRUE)
            split.data = split(x,strato)
            if(method==2)
              ans = lapply(split.data,setKgrp2)
            else
              ans = lapply(split.data,setKgrp)
            ans = do.call("rbind",ans)
            return(ans)

          }

        out = lapply(obj.split,sampleStrata)
      }
    else ## if no stratify variable is supplied
      {
        if(method==2)
          out = lapply(obj.split,setKgrp2)
        else
          out = lapply(obj.split,setKgrp)

      }

    out = do.call("rbind",out)
    object$Kgrp <- out[order(out[,1]) , 2]
    object <- object[, -which(colnames(object) %in% c("ind"))]
    return(object)
  }

