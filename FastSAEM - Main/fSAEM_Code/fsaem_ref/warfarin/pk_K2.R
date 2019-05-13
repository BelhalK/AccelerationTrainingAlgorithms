
  setwd("/Users/karimimohammedbelhal/Desktop/CSDA_code/Dir")
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_cov.R') 
  source('func_distcond.R') 
  source('func_FIM.R') 
  source('func_ggplot2.R') 
  source('func_plots.R') 
  source('func_simulations.R') 
  source('ggplot2_global.R') 
  # source('KL.R') 
  #source('vi.R') 
  source('global.R')
  source('main.R')
  source('mcmc_main.R') 
  source('main_estep.R')
  source('main_estep_mcmc.R') 
  source('main_estep_morekernels.R') 
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('plots_ggplot2.R') 
  source('saemix-package.R') 
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  source('main_new.R')
  source('main_estep_new2.R')
  source('main_new_mix.R')
  source('main_estep_mix.R')
  
setwd("/Users/karimimohammedbelhal/Desktop/CSDA_code/")
source("mixtureFunctions.R")


library("mlxR")
library(sgd)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#####################################################################################
# Theophylline

warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code/warfarin/warfarin_final.txt", header=T)
# warfa_data<-warfa_data[,1:4]


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

# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="warfarin"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE))


saemix.data<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

K1 = 200
K2 = 50
iterations = 1:(K1+K2+1)
end = K1+K2

#RWM
options<-list(seed=395246,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.burn =0)
theo_ref<-data.frame(saemix_new(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

graphConvMC_twokernels(theo_ref,theo_ref, title="new kernel")
#ref (map always)
theo_ref[end,]
theo_new_ref[end,]
options.new<-list(seed=395246,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,6),nbiter.saemix = c(K1,K2),nbiter.burn =0)
theo_new_ref<-data.frame(saemix_new(saemix.model,saemix.data,options.new))
theo_new_ref <- cbind(iterations, theo_new_ref)

graphConvMC_twokernels(theo_ref,theo_new_ref, title="new kernel")


#RWM vs always MAP (ref)



replicate = 1
seed0 = 395246

#RWM
final_rwm <- 0
final_mix <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(1,7,1,0,0,0),c(0.8,7.2,0.8,0,0,0),c(1.2,6.8,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model<-saemixModel(model=model1cpt,description="warfarin"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE))


  options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,nbiter.burn =0)
  theo_ref<-data.frame(saemix_new(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,theo_ref)

  options.new<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,6),nbiter.saemix = c(K1,K2),nbiter.burn =0)
  theo_mix<-data.frame(saemix_new(saemix.model,saemix.data,options.new))
  theo_mix <- cbind(iterations, theo_mix)
  theo_mix['individual'] <- m
  final_mix <- rbind(final_mix,theo_mix)
}

graphConvMC_new(final_rwm, title="RWM")
graphConvMC_diff(final_rwm,final_mix, title="RWM")
# graphConvMC_twokernels(final_rwm[final_rwm$individual==2,],final_rwm[final_rwm$individual==1,], title="EM")






