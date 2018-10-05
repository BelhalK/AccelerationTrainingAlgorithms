setwd("/Users/karimimohammedbelhal/Desktop/ongoing_research/mala/R")
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

setwd("/Users/karimimohammedbelhal/Desktop/ongoing_research/mala")

source('mixtureFunctions.R') 
library('rCMA')
###WARFA



K1 = 200
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2





warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/ongoing_research/mala/data/warfarin_data.txt", header=T)
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



options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_optim<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_with))
warfa_optim <- cbind(iterations, warfa_optim)

options_warfa_mala<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,0,0,0,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0, sigma.val=0.001,gamma.val=0.1)
warfa_mala<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_mala))
warfa_mala <- cbind(iterations, warfa_mala)



#RWM
final_optim <- 0
final_mala <- 0
replicate = 3
seed0 = 395246

for (m in 1:replicate){
  print(m)
  l = list(c(1,7,1,0,0,0),c(0.8,7.2,0.8,0,0,0),c(1.2,6.8,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(3,15,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

  #No var
  ##### Optim (fmin search)
  options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
  warfa_optim<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_with))
  warfa_optim <- cbind(iterations, warfa_optim)
  warfa_optim['individual'] <- m
  final_optim <- rbind(final_optim,warfa_optim)
  ##### MALA
  options_warfa_mala<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,0,0,0,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0, sigma.val=m*0.01,gamma.val=m*0.1)
  warfa_mala<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_mala))
  warfa_mala <- cbind(iterations, warfa_mala)
  warfa_mala['individual'] <- m
  final_mala <- rbind(final_mala,warfa_mala)
}

graphConvMC_diff(final_optim,final_mala, title="")
# graphConvMC_diff4(final_optim,final_optim,final_optim,final_optim, title="")

#black: optim
#blue: av
#red: av newkernel
#green: bayes


graphConvMC_diff(final_av,final_avnew, title="")

