library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
# save.image("realwarfa_mcmc_conv_varwithnew.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/VariationalInference/variationalSAEM/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('estep_mcmc.R')
  source('variationalinferencelinear.R')
  source('main.R')
  source('main_estep.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('check_linearvslaplace.R')
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  # source('SaemixRes_c.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  source('graphplot.R')

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/VariationalInference/variationalSAEM")


require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

K1 = 400
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

# Doc
oxboys.saemix<-read.table( "data/ox_synth.csv",header=T,na=".",sep=",")
oxboys.saemix_less <- oxboys.saemix[,]
n <- length(unique(oxboys.saemix_less$id))

saemix.data<-saemixData(name.data=oxboys.saemix_less,header=TRUE,
  name.group=c("id"),name.predictors=c("time"),name.response=c("y"),
  units=list(x="yr",y="cm"))

# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

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

saemix.model<-saemixModel(model=growth.linear,description="Linear model",type="structural",
  psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2,byrow=TRUE),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE), 
  error.model="constant")


L_mcmc=500
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0,0,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model,saemix.data,options_warfa)$eta_ref

options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model,saemix.data,options_warfanew)$eta

options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=2,nbiter.mcmc = c(0,0,0,6,0,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
Gamma<-mcmc(saemix.model,saemix.data,options_warfanew)$Gamma

# options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=2,nbiter.mcmc = c(0,0,0,6,0,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# test<-check.mu.gamma(saemix.model,saemix.data,options_warfanew)

K=100
variational.post.options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.gd = c(K),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),Gamma.laplace=Gamma)
variational.post<-variational.inference.linear(saemix.model,saemix.data,variational.post.options)
variational.post$mu



options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=variational.post$mu,
        Gamma = variational.post$Gamma,
        nbiter.mcmc = c(0,0,0,0,0,0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),
        nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
vi<-mcmc(saemix.model,saemix.data,options_warfavi)$eta


start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 2))

i = 10
indetabarref <- ref[[i]]
indexpecref <- data.frame(apply(indetabarref[-(1:start_interval),], 2, cummean))
indexpecref$iteration <- 1:(L_mcmc-start_interval)


indsdref <- 0
indvar <- data.frame(apply(ref[[i]][-(1:start_interval),]^2, 2, cummean))
indmeansq <- data.frame(apply(ref[[i]][-(1:start_interval),], 2, cummean))^2
indsdref <- indsdref + sqrt(pmax(zero,indvar - indmeansq))
indsdref$iteration <- 1:(L_mcmc-start_interval)


indetabarnew <- new[[i]]
indexpecnew <- data.frame(apply(indetabarnew[-(1:start_interval),], 2, cummean))
indexpecnew$iteration <- 1:(L_mcmc-start_interval)


indsdnew <- 0
indvar <- data.frame(apply(new[[i]][-(1:start_interval),]^2, 2, cummean))
indmeansq <- data.frame(apply(new[[i]][-(1:start_interval),], 2, cummean))^2
indsdnew <- indsdnew + sqrt(pmax(zero,indvar - indmeansq))
indsdnew$iteration <- 1:(L_mcmc-start_interval)

indetabarvi <- vi[[i]]
indexpecvi <- data.frame(apply(indetabarvi[-(1:start_interval),], 2, cummean))
indexpecvi$iteration <- 1:(L_mcmc-start_interval)


indsdvi <- 0
indvar <- data.frame(apply(vi[[i]][-(1:start_interval),]^2, 2, cummean))
indmeansq <- data.frame(apply(vi[[i]][-(1:start_interval),], 2, cummean))^2
indsdvi <- indsdvi + sqrt(pmax(zero,indvar - indmeansq))
indsdvi$iteration <- 1:(L_mcmc-start_interval)

plotmcmc(indexpecref[,c(3,1:2)],indexpecnew[,c(3,1:2)],title=paste("mean",i))
plotconv3(indexpecref[,c(3,1:2)],indexpecnew[,c(3,1:2)],indexpecvi[,c(3,1:2)],title="mean")