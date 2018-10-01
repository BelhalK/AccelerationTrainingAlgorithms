setwd("/Users/karimimohammedbelhal/Desktop/paramToRV/saemixnonvariability")
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
setwd("/Users/karimimohammedbelhal/Desktop/paramToRV/")
source('plots.R') 

library('rCMA')
###WARFA
warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/paramToRV/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]
  CL<-k*V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

saemix.model_warfanovar<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE))


K1 = 150
K2 = 20
iterations = 1:(K1+K2+1)
end = K1+K2


#With var no sa
options_warfa_without<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_without<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_without))
warfa_without <- cbind(iterations, warfa_without)

graphConvMC_twokernels(warfa_without,warfa_without)


options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
warfa_newkernel<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_newkernel))
warfa_newkernel <- cbind(iterations, warfa_newkernel)

graphConvMC_twokernelslog(warfa_without,warfa_newkernel)
graphConvMC_twokernels(warfa_without,warfa_newkernel)


#With var sa
options_warfa_without_with_sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_without_sa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_without_with_sa))
warfa_without_sa <- cbind(iterations, warfa_without_sa)

graphConvMC_twokernels(warfa_without,warfa_without_sa)
graphConvMC_twokernelslog(warfa_without,warfa_without_sa)


options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
warfa_newkernel_sa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_newkernel))
warfa_newkernel_sa <- cbind(iterations, warfa_newkernel_sa)

graphConvMC_twokernelslog(warfa_without_sa,warfa_newkernel_sa)


#No var no sa
options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_withnosa<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
warfa_withnosa <- cbind(iterations, warfa_withnosa)

options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_withnosanocma<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
warfa_withnosanocma <- cbind(iterations, warfa_withnosanocma)

graphConvMC_twokernels(warfa_withnosa,warfa_withnosanocma)

graphConvMC_twokernels(warfa_withnosa,warfa_without)

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
warfa_newkernelnosa<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_newkernel))
warfa_newkernelnosa <- cbind(iterations, warfa_newkernelnosa)

graphConvMC_twokernelslog(warfa_withnosa,warfa_newkernelnosa)

#No var sa
options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_withsa<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
warfa_withsa <- cbind(iterations, warfa_withsa)

graphConvMC_twokernels(warfa_withav,warfa_without)

graphConvMC_twokernelslog(warfa_newkernel_sa[,1:6],warfa_newkernelsa[,1:6])

graphConvMC_twokernelslog(warfa_without_sa[,1:6],warfa_withsa[,1:6])

graphConvMC_twokernelslog(warfa_without_sa,warfa_without_sa)
graphConvMC_twokernelslog(warfa_withsa,warfa_withsa)

graphConvMC_twokernelslog(warfa_newkernel_sa,warfa_newkernel_sa)
graphConvMC_twokernelslog(warfa_newkernelsa,warfa_newkernelsa)

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
warfa_newkernelsa<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_newkernel))
warfa_newkernelsa <- cbind(iterations, warfa_newkernelsa)

graphConvMC_twokernelslog(warfa_withsa,warfa_newkernelsa)

#No var av
options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=1)
warfa_withav<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
warfa_withav <- cbind(iterations, warfa_withav)

graphConvMC_twokernels(warfa_withav,warfa_without)


options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=1)
warfa_newkernelav<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_newkernel))
warfa_newkernelav <- cbind(iterations, warfa_newkernelav)

graphConvMC_twokernelslog(warfa_withav,warfa_newkernelav)

#AV=1
graphConvMC_twokernelslog(warfa_withav,warfa_newkernelav, title="AV")
#Simulate annealing
graphConvMC_twokernelslog(warfa_without_sa,warfa_newkernel_sa, title="SA")
