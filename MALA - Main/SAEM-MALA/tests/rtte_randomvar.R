setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/saemixrandomvariable")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('main.R')
  source('main_estep.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/")
source('plots.R') 

###WARFA
###RTTE
timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/data/rtte_data.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
saemix.data_rtte<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
timetoevent.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2]
  N <- nrow(psi)
  Nj <- length(T)
  censoringtime = 20
  lambda <- psi[id,1]
  beta <- psi[id,2]
  init <- which(T==0)
  cens <- which(T==censoringtime)
  ind <- setdiff(1:Nj, append(init,cens))
  hazard <- (beta/lambda)*(T/lambda)^(beta-1)
  H <- (T/lambda)^beta
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens] + H[cens-1]
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
  return(logpdf)
}

saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

saemix.model_rttenovar<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, 
  byrow=TRUE))

##RUNS

K1 = 150
K2 = 20
iterations = 1:(K1+K2+1)
end = K1+K2



#with var no sa
options_rtte<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2), nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, av=0)
rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
rtte <- cbind(iterations, rtte)


options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.saemix = c(K1,K2), nbiter.sa=0,displayProgress=TRUE,map.range=c(1:5),nbiter.burn =0, av=0)
rttenewnosa<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenew))
rttenewnosa <- cbind(iterations, rttenewnosa)

graphConvMC_twokernels(rtte,rttenewnosa)

#No var no sa but randomvariable
options_rttenosa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2), nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, av=0)
rttenosa<-data.frame(saemix(saemix.model_rttenovar,saemix.data_rtte,options_rttenosa))
rttenosa <- cbind(iterations, rttenosa)
