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
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mamyula")
source('mala_main.R')
source('main_estep_mala.R')
source("mixtureFunctions.R")

library("mlxR")
library("psych")
library("coda")
library("Matrix")

require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
iter_mcmc = 50

# Doc
theo.saemix<-read.table( "data/theo.saemix.tab",header=T,na=".")
l <- c(4.02,4.4,4.53,4.4,5.86,4,4.95,4.53,3.1,5.5,4.92,5.3)
for (i in 1:12){
  theo.saemix[(i*10-9):(i*10),'Dose'] = l[i]
}
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
# Default model, no covariate
# saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))
saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(10,10,1.05),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))

saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0,0,0))
saemix.options_mala<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc,0,0),sigma.val = 0.01,gamma.val = 0.01)
saemix.options_mamyula<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,0,iter_mcmc,0),sigma.val = 0.1,gamma.val = 0.01)
saemix.options_mamyula_nest<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,0,0,iter_mcmc),sigma.val = 0.1,gamma.val = 0.01,memory=0.001)


post_rwm<-saemix_mala(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
post_mala<-saemix_mala(saemix.model,saemix.data,saemix.options_mala)$post_mala
post_mamyula<-saemix_mala(saemix.model,saemix.data,saemix.options_mamyula)$post_mala
post_mamyula_nest<-saemix_mala(saemix.model,saemix.data,saemix.options_mamyula_nest)$post_mala


index = 4
graphConvMC_threekernels(post_rwm[[index]],post_mala[[index]],post_mamyula[[index]], title="EM")


graphConvMC_twokernels(post_rwm[[index]],post_mala[[index]], title="EM")
graphConvMC_twokernels(post_rwm[[index]],post_mamyula[[index]], title="EM")
graphConvMC_threekernels(post_rwm[[index]],post_mamyula[[index]],post_mamyula_nest[[index]], title="EM")



names(post_rwm[[index]])[2]<-paste("ka")
names(post_rwm[[index]])[3]<-paste("V")
names(post_rwm[[index]])[4]<-paste("Cl")

names(post_mala[[index]])[2]<-paste("ka")
names(post_mala[[index]])[3]<-paste("V")
names(post_mala[[index]])[4]<-paste("Cl")

graphConvMC_twokernels(post_rwm[[index]][,-c(3,4)],post_mala[[index]][,-c(3,4)])


graphConvMC_twokernels(post_rwm[[index]],post_mala[[index]], title="RWM vs MALA")
graphConvMC_twokernels(post_mala[[index]],post_mamyula[[index]], title="MALA vs mamyula")



final_rwm <- post_rwm[[1]]
for (i in 2:length(post_rwm)) {
  final_rwm <- rbind(final_rwm, post_rwm[[i]])
}


final_mala <- post_mala[[1]]
for (i in 2:length(post_mala)) {
  final_mala <- rbind(final_mala, post_mala[[i]])
}

final_mamyula <- post_mamyula[[1]]
for (i in 2:length(post_mamyula)) {
  final_mamyula <- rbind(final_mamyula, post_mamyula[[i]])
}



graphConvMC_threekernels(final_rwm,final_mala,final_mamyula, title="EM")



#Autocorrelation
rwm.obj <- as.mcmc(post_rwm[[index]])
corr_rwm <- autocorr(rwm.obj[,3])
autocorr.plot(rwm.obj[,3])

mala.obj <- as.mcmc(post_mala[[index]])
corr_mala <- autocorr(mala.obj[,3])
autocorr.plot(mala.obj[,3])

mamyula.obj <- as.mcmc(post_mamyula[[index]])
corr_mamyula <- autocorr(mamyula.obj[,3])
autocorr.plot(mamyula.obj[,3])



#MSJD
mssd(post_rwm[[index]][,3])
mssd(post_mala[[index]][,3])
mssd(post_mamyula[[index]][,3])


