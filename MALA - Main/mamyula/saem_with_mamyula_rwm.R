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
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem")

source("mixtureFunctions.R")

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mamyula")
source('mala_main.R')
source('main_estep_mala.R')
source('main_mamyula.R')
# source("mixtureFunctions.R")

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

K1 = 100
K2 = 50
iterations = 1:(K1+K2+1)
gd_step = 0.01


#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,2,0,0,0,0), nbiter.saemix = c(K1,K2))
theo_ref<-data.frame(saemix(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

#saem with mala
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,2,0,0),nbiter.saemix = c(K1,K2),sigma.val = 0.01,gamma.val=0.01)
theo_mala<-data.frame(saemix_mamyula(saemix.model,saemix.data,options.mala))
theo_mala <- cbind(iterations, theo_mala)


#saem with mamyula
options.mamyula<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,0,2,0),nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01)
theo_mamyula<-data.frame(saemix_mamyula(saemix.model,saemix.data,options.mamyula))
theo_mamyula <- cbind(iterations, theo_mamyula)

graphConvMC_twokernels(theo_ref,theo_mala, title="new kernel")
graphConvMC_twokernels(theo_ref,theo_mamyula[,-1], title="new kernel")
graphConvMC_threekernels(theo_ref,theo_mala,theo_mamyula, title="new kernel")


