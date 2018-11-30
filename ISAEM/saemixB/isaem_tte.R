
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/R")
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
  source('mixtureFunctions.R')

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/")
source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/plots.R")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)


timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/rtte3se.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
saemix.data<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))

timetoevent.model<-function(psi,id,xidep) {
T<-xidep[,1]
y<-xidep[,2]
N <- nrow(psi)
Nj <- length(T)
censoringtime = 6
lambda <- psi[id,1]
beta <- psi[id,2]
init <- which(T==0)
cens <- which(T==censoringtime)
ind <- setdiff(1:Nj, append(init,cens))
hazard <- (beta/lambda)*(T/lambda)^(beta-1)
H <- (T/lambda)^beta
logpdf <- rep(0,Nj)
logpdf[cens] <- -H[cens] + H[cens-1]
logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
return(logpdf)
}

saemix.model<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(1,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

##RUNS

K1 = 200
K2 = 100

iterations = 0:(K1+K2-1)
end = K1+K2
seed0 = 39546
map_range = c(1:10)

options<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,2),map.range=map_range,
  nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, nb.replacement=100,sampling='seq')
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])



options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,2),map.range=map_range, 
                          nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,2),map.range=map_range, nbiter.saemix = c(K1,K2),displayProgress=FALSE,
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])



theo_ref_scaled <- theo_ref
theo_mix50_scaled <- theo_mix50
theo_mix25_scaled <- theo_mix25


theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25

graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix50_scaled,theo_mix50_scaled)
# graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix85_scaled,theo_mix85_scaled)
