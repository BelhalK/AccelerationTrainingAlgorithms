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
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mala_nest")
source('main_estep_mala.R')
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mamyula")
source('main_mala.R')

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



oxboys.saemix<-read.table( "data/oxboys.saemix.tab",header=T,na=".")
oxboys.saemix_less <- oxboys.saemix
saemix.data<-saemixData(name.data=oxboys.saemix_less,header=TRUE,
  name.group=c("Subject"),name.predictors=c("age"),name.response=c("height"),
  units=list(x="yr",y="cm"))





growth.linear<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (2 columns, base and slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  base<-psi[id,1]
  slope<-psi[id,2]
  f<-base+slope*x
  return(f)
}
saemix.model<-saemixModel(model=growth.linear,description="Linear model",
  psi0=matrix(c(1,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(1,0),covariance.model=matrix(c(1,1,1,1),ncol=2,byrow=TRUE), 
  error.model="constant")


indiv = 1
seed0 = 35644
replicate = 5
iter_mcmc = 10000
burn = 400


K1 = 300
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2
gd_step = 0.00001
seed0 = 39546


#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,0,0), nbiter.saemix = c(K1,K2))
theo_ref<-data.frame(saemix_mala(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

#saem with mala
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5,0,0,0),nbiter.saemix = c(K1,K2),sigma.val = 0.01,gamma.val=0.01)
theo_mala<-data.frame(saemix_mala(saemix.model,saemix.data,options.mala))
theo_mala <- cbind(iterations, theo_mala)

graphConvMC_twokernels(theo_ref,theo_mala, title="RWM vs Laplace SAEM")

