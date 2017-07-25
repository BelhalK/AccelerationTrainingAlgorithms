#library(rstan)
setwd("/Users/karimimohammedbelhal/Desktop/variationalBayes/mcmc_R_isolate/Dir2")
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
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat")
source('laplace_main.R')
source('main_estep_laplace.R')
source("mixtureFunctions.R")

library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
iter_mcmc = 200


warfarin.saemix<-read.table("data/pkcat_data2.txt",header=T,na=".")
saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.predictors=c("amt"),name.response=c("dv","dvid"),name.covariates=c("wt","sex","age"), name.X="time")


warfarin<-function(psi,id,xidep) {
dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))


  x<-xidep[,1]
  ymax<-psi[id,1]
  xmax<-psi[id,2]
  slope<-psi[id,3]
  f<-ymax+slope*(x-xmax)
#  cat(length(f),"  ",length(ymax),"\n")
  f[x>xmax]<-ymax[x>xmax]
  return(f)
}
saemix.model<-saemixModel(model=yield.LP,description="Linear plus plateau model",   
  psi0=matrix(c(8,100,0.2,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("Ymax","Xmax","slope"))),covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")

saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0,0,0,0))
saemix.laplace<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc,0,0,0))
# saemix.fo<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,iter_mcmc,0,0))
saemix.fo2<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,0,iter_mcmc,0))
saemix.foce<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,0,0,iter_mcmc))


post_rwm<-saemix_laplace(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
post_foce<-saemix_laplace(saemix.model,saemix.data,saemix.foce)$post_newkernel
post_laplace<-saemix_laplace(saemix.model,saemix.data,saemix.laplace)$post_newkernel
# post_fo<-saemix_laplace(saemix.model,saemix.data,saemix.fo)$post_newkernel
post_fo2<-saemix_laplace(saemix.model,saemix.data,saemix.fo2)$post_newkernel


index = 1
graphConvMC_twokernels(post_rwm[[index]],post_foce[[index]], title="rwm vs foce")
# graphConvMC_twokernels(post_rwm[[index]],post_fo[[index]], title="rwm vs fo")
graphConvMC_twokernels(post_rwm[[index]],post_fo2[[index]], title="rwm vs fo2")
graphConvMC_twokernels(post_rwm[[index]],post_laplace[[index]], title="rwm vs laplace")
graphConvMC_threekernels(post_rwm[[index]],post_foce[[index]],post_laplace[[index]], title="rwm vs foce vs laplace")
graphConvMC_threekernels(post_rwm[[index]],post_foce[[index]],post_fo2[[index]], title="rwm vs foce vs laplace")


final_rwm <- post_rwm[[1]]
for (i in 2:length(post_rwm)) {
  final_rwm <- rbind(final_rwm, post_rwm[[i]])
}


final_foce <- post_foce[[1]]
for (i in 2:length(post_foce)) {
  final_foce <- rbind(final_foce, post_foce[[i]])
}


final_laplace <- post_laplace[[1]]
for (i in 2:length(post_laplace)) {
  final_laplace <- rbind(final_laplace, post_laplace[[i]])
}

final_fo2 <- post_fo2[[1]]
for (i in 2:length(post_fo2)) {
  final_fo2 <- rbind(final_fo2, post_fo2[[i]])
}
graphConvMC_new(final_fo2, title="VB Linear case")

#ALl individual posteriors
graphConvMC_new(final_rwm, title="RWM")
graphConvMC_new(final_laplace, title="VB Linear case")

#first individual posteriors
graphConvMC_new(post_rwm[[index]], title="EM")

graphConvMC_twokernels(final_rwm,final_foce, title="EM")
graphConvMC_twokernels(final_rwm,final_laplace, title="EM")
graphConvMC_threekernels(final_rwm,final_fo2,final_foce, title="EM")





#Autocorrelation
rwm.obj <- as.mcmc(post_rwm[[1]])
corr_rwm <- autocorr(rwm.obj[,2])
autocorr.plot(rwm.obj[,2])

foce.obj <- as.mcmc(post_foce[[1]])
corr_foce <- autocorr(foce.obj[,2])
autocorr.plot(foce.obj[,2])

fo2.obj <- as.mcmc(post_fo2[[1]])
corr_fo2 <- autocorr(fo2.obj[,2])
autocorr.plot(fo2.obj[,2])

laplace.obj <- as.mcmc(post_laplace[[1]])
corr_laplace <- autocorr(laplace.obj[,2])
autocorr.plot(laplace.obj[,2])

#MSJD
mssd(post_rwm[[index]][,2])
mssd(post_foce[[index]][,2])
mssd(post_fo2[[index]][,2])
mssd(post_laplace[[index]][,2])



