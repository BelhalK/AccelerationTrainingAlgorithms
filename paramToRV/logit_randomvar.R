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

library('rCMA')
###logit


logit_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/data/Data_logistic_1D_simule.txt", header=T)
saemix.data_logit<-saemixData(name.data=logit_data,header=TRUE,sep=" ",na=NA, name.group=c("group"),
  name.predictors=c("X"),name.response=c("Y"), name.X="X")

model1cpt<-function(psi,id,xidep) { 
  tim<-xidep[,1]  
  p0<-psi[id,1]
  alpha<-psi[id,2]
  tau<-psi[id,3]
  ypred<-1/(1+((1/p0)-1)*exp(-alpha*(tim-tau)/(p0*(1-p0))))
  return(ypred)
}


saemix.model_logit<-saemixModel(model=model1cpt,description="logitrin",type="structural"
  ,psi0=matrix(c(1,7,1),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

saemix.model_logitnovar<-saemixModel(model=model1cpt,description="logitrin",type="structural"
  ,psi0=matrix(c(1,7,1),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE))


K1 = 500
K2 = 500
iterations = 1:(K1+K2+1)
end = K1+K2


#With var no sa
options_logit_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
logit_without<-data.frame(saemix(saemix.model_logit,saemix.data_logit,options_logit_without))
logit_without <- cbind(iterations, logit_without)

graphConvMC_twokernels(logit_without,logit_without)


options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
logit_newkernel<-data.frame(saemix(saemix.model_logit,saemix.data_logit,options_newkernel))
logit_newkernel <- cbind(iterations, logit_newkernel)


#no var no sa
options_logit_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
logit_without<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_logit_without))
logit_without <- cbind(iterations, logit_without)

graphConvMC_twokernels(logit_without,logit_without)


options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
logit_newkernel<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_newkernel))
logit_newkernel <- cbind(iterations, logit_newkernel)


#No var no sa but randomvariable
options_logit_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0,map.range=c(0))
logit_withnosa<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_logit_with))
logit_withnosa <- cbind(iterations, logit_withnosa)
graphConvMC_twokernels(logit_without[,1:6],logit_withnosa)

options_logit_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0,map.range=c(0))
logit_withnosa<-data.frame(saemix(saemix.model_logitnovar_init,saemix.data_logit,options_logit_with))
logit_withnosa <- cbind(iterations, logit_withnosa)