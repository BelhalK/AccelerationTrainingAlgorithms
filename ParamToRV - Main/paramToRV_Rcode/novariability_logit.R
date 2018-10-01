setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/saemixnonvariability")
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
  ,psi0=matrix(c(2,2,2),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

saemix.model_logitnovar<-saemixModel(model=model1cpt,description="logitrin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE))


K1 = 300
K2 = 300
iterations = 1:(K1+K2+1)
end = K1+K2


#With var no sa
options_logit_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
logit_without<-data.frame(saemix(saemix.model_logit,saemix.data_logit,options_logit_without))
logit_without <- cbind(iterations, logit_without)


graphConvMC_twokernels(logit_without,logit_without)


options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:20), av=0)
logit_newkernel<-data.frame(saemix(saemix.model_logit,saemix.data_logit,options_newkernel))
logit_newkernel <- cbind(iterations, logit_newkernel)

graphConvMC_twokernelslog(logit_without,logit_newkernel)
graphConvMC_twokernels(logit_without,logit_newkernel)


#No var no sa
options_logit_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
logit_withnosa<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_logit_with))
logit_withnosa <- cbind(iterations, logit_withnosa)

options_logit_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
logit_withnosanocma<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_logit_with))
logit_withnosanocma <- cbind(iterations, logit_withnosanocma)

graphConvMC_twokernels(logit_withnosa,logit_withnosanocma)

graphConvMC_twokernels(logit_withnosa,logit_without)

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
logit_newkernelnosa<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_newkernel))
logit_newkernelnosa <- cbind(iterations, logit_newkernelnosa)

graphConvMC_twokernelslog(logit_withnosa,logit_newkernelnosa)

#No var sa
options_logit_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
logit_withsa<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_logit_with))
logit_withsa <- cbind(iterations, logit_withsa)

graphConvMC_twokernels(logit_withav,logit_without)

graphConvMC_twokernelslog(logit_newkernel_sa[,1:6],logit_newkernelsa[,1:6])

graphConvMC_twokernelslog(logit_without_sa[,1:6],logit_withsa[,1:6])

graphConvMC_twokernelslog(logit_without_sa,logit_without_sa)
graphConvMC_twokernelslog(logit_withsa,logit_withsa)

graphConvMC_twokernelslog(logit_newkernel_sa,logit_newkernel_sa)
graphConvMC_twokernelslog(logit_newkernelsa,logit_newkernelsa)

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
logit_newkernelsa<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_newkernel))
logit_newkernelsa <- cbind(iterations, logit_newkernelsa)

graphConvMC_twokernelslog(logit_withsa,logit_newkernelsa)

#No var av
options_logit_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=1)
logit_withav<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_logit_with))
logit_withav <- cbind(iterations, logit_withav)

graphConvMC_twokernels(logit_withav,logit_without)


options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:3), av=1)
logit_newkernelav<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_newkernel))
logit_newkernelav <- cbind(iterations, logit_newkernelav)

graphConvMC_twokernelslog(logit_withav,logit_newkernelav)

#AV=1
graphConvMC_twokernelslog(logit_withav,logit_newkernelav, title="AV")
#Simulate annealing
graphConvMC_twokernelslog(logit_without_sa,logit_newkernel_sa, title="SA")
