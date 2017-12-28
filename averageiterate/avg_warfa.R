setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/averageiterate")
source('mixtureFunctions.R') 

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/averageiterate/avg")
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
  source('main_initialiseMainAlgoavg.R')
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 


###WARFA
warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/averageiterate/data/warfarin_data.txt", header=T)
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



K1 = 1000
K2 = 3000
iterations = 1:(K1+K2+1)
end = K1+K2


#With var no sa
options.ref<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=0)
warfa.ref<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options.ref))
warfa.ref <- cbind(iterations, warfa.ref)
# graphConvMC_twokernels(warfa.ref,warfa.ref)

options.avg<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=1)
warfa.avg<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options.avg))
warfa.avg <- cbind(iterations, warfa.avg)

# graphConvMC_twokernels(warfa.avg,warfa.avg)
graphConvMC_twokernels(warfa.ref,warfa.avg)

graphConvMC_twokernels(warfa.ref[K1:end,],warfa.avg[K1:end,])

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(1:5), av=0,avg=0)
warfa_newkernel<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_newkernel))
warfa_newkernel <- cbind(iterations, warfa_newkernel)

options_newavg<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(1:5), av=0,avg=1)
warfa_newavg<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_newavg))
warfa_newavg <- cbind(iterations, warfa_newavg)
graphConvMC_twokernels(warfa_without,warfa_newkernel)




