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
  # source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('plots_ggplot2.R') 
  source('saemix-package.R') 
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  # source('SaemixObject.R') 
  source('zzz.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat")
source('main_cat2.R')
source('main_estep_cat2.R')
source('main_estep_cat.R')
source('main_mstep_cat.R') 
source('post_cat.R') 
source('func_aux_cat.R') 
source('SaemixObject_cat.R') 
source('main_initialiseMainAlgo_cat.R') 
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


# cat_data.saemix<-read.table("data/categorical1_data.txt",header=T,na=".")
# saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("ID"),name.response=c("Y"),name.predictors=c("Y"), name.X=c("TIME"))



# cat_data.saemix<-read.table("data/categorical1_data_less.txt",header=T,na=".")
# cat_data.saemix<-read.table("data/categorical1_data_less2.txt",header=T,na=".")
# cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat1.csv", header=T, sep=",")
cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat2.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("y"), name.X=c("time"))



cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]

th1 <- psi[id,1]
th2 <- psi[id,2]
th3 <- psi[id,3]

P0 <- 1/(1+exp(-th1))
Pcum1 <- 1/(1+exp(-th1-th2))
Pcum2 <- 1/(1+exp(-th1-th2-th3))


P1 <- Pcum1 - P0
P2 <- Pcum2 - Pcum1
P3 <- 1 - Pcum2

P.obs = (level==0)*P0+(level==1)*P1+(level==2)*P2+(level==3)*P3

return(P.obs)
}


saemix.model<-saemixModel(model=cat_data.model,description="cat model",   
  psi0=matrix(c(1,1,1),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))),covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3, 
  byrow=TRUE),error.model="constant")



gd_step = 0.01
end = K1+K2
seed0 = 444
iter_mcmc = 2000

#RWM
theo_ref <- NULL
saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0))
theo_ref<-saemix_post_cat(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
indiv = 1
matplot ((1:iter_mcmc), cbind (theo_ref[[indiv]][,1], theo_ref[[indiv]][,1]), pch = 19,type="l", col=c(2,4))
legend("bottomleft", inset=.05, legend=c("ref", "new"), pch=1, col=c(2,4), horiz=TRUE)

#MAP then RWM
cat_saem <- NULL
saemix.foce<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc),map.range=c(1))
cat_saem<-saemix_post_cat(saemix.model,saemix.data,saemix.foce)$post_new


indiv = 3
matplot ((1:iter_mcmc), cbind (theo_ref[[indiv]][,1], cat_saem[[indiv]][,1]), pch = 19,type="l", col=c(2,4))
legend("bottomleft", inset=.05, legend=c("ref", "new"), pch=1, col=c(2,4), horiz=TRUE)
