genData <- function(n,p,intervene=NULL) {
  
  # exogenous variables
  for(i in 1:2){
    assign(paste0("U.Lt",i), rnorm(n, mean=0, sd=1))
    assign(paste0("U.At",i), rnorm(n,0,1))
    assign(paste0("U.Yt",i+1), rnorm(n,0,1))

  }
  
  # endogenous variables
  L1 = as.numeric(U.Lt1 < 0.4) #most people controlled at BL
  if(is.null(intervene)){
    A1 = as.numeric( U.At1+1 < plogis(L1))#people controlled/healthier more likely to go on GLP1
  }else{
    A1<-intervene
  }
  Y2 = as.numeric(U.Yt2 > plogis((L1+(A1*2))))#less people with tx+controlled L1 have outcome
  L2 = ifelse(Y2==0, as.numeric(plogis(A1 + L1 + U.Lt2) > 0.4),0)#NA
  if(is.null(intervene)){
    A2 = ifelse(Y2==0,as.numeric(U.At2+1 < plogis((L2)+A1+L1)),A1)#NA
  }else{
    A2<-intervene
  }
  Y3 = ifelse(Y2==0, as.numeric(U.Yt3+1 < plogis(L1+A1+(L2)+A2)),1)#NA
  # dataframe with endogenous variables
  ObsData <- data.table(L1,A1,Y2,L2,A2,Y3)
  table(is.na(ObsData$L2))

  return(ObsData)
}

run_analysis <- function(data,Lnodes,Cnodes,Ynodes=c("Y2","Y3"),str_name,deterministic.g.function=NULL,abar=list(c(1,1),c(0,0))){
  
  analysis <- ltmle(data,
                    Anodes=c("A1","A2"),
                    Lnodes=Lnodes,
                    Cnodes=Cnodes,
                    Ynodes=Ynodes,
                    abar= abar,
                    survivalOutcome=TRUE,
                    deterministic.g.function=deterministic.g.function)
  res <-summary(analysis)
  res.ate <- as.data.frame(c(res$effect.measures$ATE,str_name)) %>%
    rename(ate.long.name=long.name,
           ate=estimate,
           ate.sd=std.dev,
           ate.pval=pvalue,
           ate.ci.lb=CI.2.5.,
           ate.ci.ub=  CI.97.5.,
           ate.log.std.err=log.std.err,
           analysis=(paste0("X.",str_name,".")))
  return(res.ate)
}
