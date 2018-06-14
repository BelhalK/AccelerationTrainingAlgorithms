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


setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/pd")

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

library(saemix)
PD1.saemix<-read.table( "PD1.saemix.tab",header=T,na=".")
PD2.saemix<-read.table( "PD2.saemix.tab",header=T,na=".")
saemix.data1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))
saemix.data2<-saemixData(name.data=PD2.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))
modelemax<-function(psi,id,xidep) {
# input:
# psi : matrix of parameters (3 columns, E0, Emax, EC50)
# id : vector of indices
# xidep : dependent variables (same nb of rows as length of id)
# returns:
# a vector of predictions of length equal to length of id
dose<-xidep[,1]
e0<-psi[id,1]
emax<-psi[id,2]
e50<-psi[id,3]
f<-e0+emax*dose/(e50+dose)
return(f)
}
saemix.model<-saemixModel(model=modelemax,description="Emax model",
psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")


K1 = 100
K2 = 50
iterations = 1:(K1+K2+1)
gd_step = 0.01


#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,0), nbiter.saemix = c(K1,K2))
theo_ref<-data.frame(saemix(saemix.model,saemix.data2,options))
theo_ref <- cbind(iterations, theo_ref)

#saem with mala
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5,0,0),nbiter.saemix = c(K1,K2),sigma.val = 0.01,gamma.val=0.01)
theo_mala<-data.frame(saemix_mamyula(saemix.model,saemix.data2,options.mala))
theo_mala <- cbind(iterations, theo_mala)


#saem with mamyula
options.mamyula<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,0,5,0),nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01)
theo_mamyula<-data.frame(saemix_mamyula(saemix.model,saemix.data2,options.mamyula))
theo_mamyula <- cbind(iterations, theo_mamyula)

graphConvMC_twokernels(theo_ref,theo_mala, title="new kernel")
graphConvMC_twokernels(theo_ref,theo_mamyula, title="new kernel")
graphConvMC_threekernels(theo_ref,theo_mala,theo_mamyula, title="new kernel")


