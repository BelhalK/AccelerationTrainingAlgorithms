library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Desktop/csda_new/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('estep_mcmc.R')
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
  
setwd("/Users/karimimohammedbelhal/Desktop/csda_new")
source('graphplot.R') 
warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/csda_new/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")
n <- length(unique(warfa_data$id))
model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]

  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


##RUNS

K1 = 400
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
warfa <- cbind(iterations, warfa)


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(K1))
warfanew<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew))
warfanew <- cbind(iterations, warfanew)


graphConvMC_twokernels(warfa,warfanew)


#compareMCMC



# saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(0.611,7.61,0.0178),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(0.43,0,0,0,0.04,0,0,0,0.06),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.5,7,0.1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


L_mcmc=2000
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta_ref

options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta


start_interval <- 100
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 3))

etabarref <- 1/n*Reduce("+",ref)
expecref <- data.frame(apply(etabarref[-(1:start_interval),], 2, cummean))
expecref$iteration <- 1:(L_mcmc-start_interval)


sdref <- 0
for (i in 1:n){
  var <- data.frame(apply(ref[[i]][-(1:start_interval),]^2, 2, cummean))
  meansq <- data.frame(apply(ref[[i]][-(1:start_interval),], 2, cummean))^2
  sdref <- sdref + sqrt(pmax(zero,var - meansq))
}


sdref <- 1/n*sdref
sdref$iteration <- 1:(L_mcmc-start_interval)


etabarnew <- 1/n*Reduce("+",new)
expecnew <- data.frame(apply(etabarnew[-(1:start_interval),], 2, cummean))
expecnew$iteration <- 1:(L_mcmc-start_interval)


sdnew <- 0
for (i in 1:n){
  var <- data.frame(apply(new[[i]][-(1:start_interval),]^2, 2, cummean))
  meansq <- data.frame(apply(new[[i]][-(1:start_interval),], 2, cummean))^2
  sdnew <- sdnew + sqrt(pmax(zero,var - meansq))
}

sdnew <- 1/n*sdnew
sdnew$iteration <- 1:(L_mcmc-start_interval)


plotmcmc(expecref[,c(4,1:3)],expecnew[,c(4,1:3)],title="mean")
plotmcmc(sdref[,c(4,1:3)],sdnew[,c(4,1:3)],title="sd")


etaref <- 1/n*Reduce("+",ref)
etaref$iteration <- 1:(L_mcmc)
# plotmcmc(etaref[,c(4,1:3)],etaref[,c(4,1:3)],title="mean")

etanew <- 1/n*Reduce("+",new)
etanew$iteration <- 1:(L_mcmc)

plotmcmc(etaref[,c(4,1:3)],etanew[,c(4,1:3)],title="mean")



for (i in 6:12){
ref[[i]]$iteration <- 1:(L_mcmc)
new[[i]]$iteration <- 1:(L_mcmc)
plotmcmc(ref[[i]][,c(4,1:3)],new[[i]][,c(4,1:3)],title="mean")
}

plotmcmc(ref[[9]][,c(4,1:3)],new[[9]][,c(4,1:3)],title="mean")
plotmcmc(ref[[5]][,c(4,1:3)],new[[5]][,c(4,1:3)],title="mean")

# options_warfamix<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(iter_mcmc,iter_mcmc,0,iter_mcmc),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# mix<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfamix)$eta

# meanmix <- 1/n*Reduce("+",mix)
# meanmix$iteration <- 1:(3*iter_mcmc)


# expecmix <- colMeans(meanmix[-(1:100),2:4])
# sdmix <- meanmix
# sdmix[,2:4] <- sqrt((meanmix[,2:4]-expecmix)*(meanmix[,2:4]-expecmix))



# plotmcmc(meanmix[,1:4],meanref[,1:4],title="mean")
# plotmcmc(sdmix[,1:4],sdref[,1:4],title="sd")



