# load("warfa_isaem_newkernel.RData")
# save.image("warfa_isaem_newkernel.RData")
# save.image("warfa_isaem_LARGEN.RData")

####20 CHAINS
load("warfa_isaem.RData")
# save.image("warfa_isaem.RData")

# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incremental/R")
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
  # source('SaemixRes_c.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  
  source('main_incremental.R')
  source('main_estep_incremental.R')


  source('/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incremental/R/mixtureFunctions.R')
  source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incremental/plots.R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incremental")


library("mlxR")
library("rlist")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)


# warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/csda_new/data/warfarin_data.txt", header=T)



model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  k <- Cl/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*time)-exp(-ka*time))
  return(ypred)
}

warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incremental/data/warfarin_data.txt", header=T)
saemix.data<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")


saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(3,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))



K1 = 200
K2 = 30
iterations = 0:(K1+K2-1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50

seed0=3456

nchains = 1
gamma = 1
options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains,nbiter.mcmc = c(2,2,2,0), 
  nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, 
  map.range=c(0), nb.replacement=100,sampling='randomiter',gamma=gamma)
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])


# options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, 
#   nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),
#   nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randomiter',gamma=gamma)
# theo_mix75<-saemix_incremental(saemix.model,saemix.data,options.incremental75)
# theo_mix75 <- data.frame(theo_mix75$param)
# theo_mix75 <- cbind(iterations, theo_mix75[-1,])


# options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,0), 
#                           nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,
#                           nbiter.burn =0, nb.replacement=50,sampling='seq',gamma=gamma)
# theo_mix50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
# theo_mix50 <- data.frame(theo_mix50$param)
# theo_mix50 <- cbind(iterations, theo_mix50[-1,])


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=25,sampling='randompass',gamma=gamma)
theo_mix25online<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25online <- data.frame(theo_mix25online$param)
theo_mix25online <- cbind(iterations, theo_mix25online[-1,])

theo_mix25online_scaled <- theo_mix25online
theo_mix25online_scaled$iterations = theo_mix25online_scaled$iterations*0.25


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='seq',gamma=gamma)
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])


theo_ref_scaled <- theo_ref
theo_mix25_scaled <- theo_mix25
theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25



graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix25_scaled)
graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix25online_scaled)
