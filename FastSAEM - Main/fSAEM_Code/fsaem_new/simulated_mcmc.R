library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
load("sim_mcmc_conv.RData")
# save.image("sim_mcmc_conv.RData")
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


ka_true <- 1
V_true <- 8
k_true <- 0.01
o_ka <- 0.5
o_V <- 0.2
o_k <- 0.2
a_true <- 1

myModel <- inlineModel("


[INDIVIDUAL]
input = {ka_pop, V_pop, k_pop, omega_ka, omega_V, omega_k}
DEFINITION:
ka = {distribution=lognormal, reference=ka_pop, sd=omega_ka}
V  = {distribution=lognormal, reference=V_pop,  sd=omega_V }
k = {distribution=lognormal, reference=k_pop, sd=omega_k}


[LONGITUDINAL]
input = {ka, V, k,a}
EQUATION:
C = pkmodel(ka,V,k)
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=50

pop.param   <- c(
  ka_pop  = ka_true,    omega_ka  = o_ka,
  V_pop   = V_true,   omega_V   = o_V,
  k_pop  = k_true,    omega_k  = o_k, a =a_true)
  
res <- simulx(model     = myModel,
              parameter = pop.param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(0,20,by=1)))
  

  warfarin.saemix <- res$y
  warfarin.saemix$amount <- 100
writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Desktop/mcmc_pk.csv")


saemix.data_warfa<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")


n <- length(unique(warfarin.saemix$id))
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



saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.8,7.31,0.0278),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(0.2,0,0,0,0.18,0,0,0,0.03),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

# saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(0.5,7,0.1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))


L_mcmc=50000
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta_ref

options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta


start_interval <- 200
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
plotmcmc(sdref[-c(1:10),c(4,1:3)],sdnew[-c(1:10),c(4,1:3)],title="sd")


etaref <- 1/n*Reduce("+",ref)
etaref$iteration <- 1:(L_mcmc)
# plotmcmc(etaref[,c(4,1:3)],etaref[,c(4,1:3)],title="mean")

etanew <- 1/n*Reduce("+",new)
etanew$iteration <- 1:(L_mcmc)

plotmcmc(etaref[,c(4,1:3)],etanew[,c(4,1:3)],title="mean")



for (i in 1:5){
ref[[i]]$iteration <- 1:(L_mcmc)
new[[i]]$iteration <- 1:(L_mcmc)
plotmcmc(ref[[i]][,c(4,1:3)],new[[i]][,c(4,1:3)],title="mean")
}

plotmcmc(ref[[9]][,c(4,1:3)],new[[9]][,c(4,1:3)],title="mean")
plotmcmc(ref[[5]][,c(4,1:3)],new[[5]][,c(4,1:3)],title="mean")

#one invdiv


start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 3))


for (i in 1:2){
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


plotmcmc(indexpecref[,c(4,1:3)],indexpecnew[,c(4,1:3)],title=paste("mean",i))
plotmcmc(indsdref[-c(1:10),c(4,1:3)],indsdnew[-c(1:10),c(4,1:3)],title=paste("sd",i))
}


start_q <- 1

i <- 5
qref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[[i]][1:k,dim], 0.05)
    qref[[dim]][k,2] <- quantile(ref[[i]][1:k,dim], 0.5)
    qref[[dim]][k,3] <- quantile(ref[[i]][1:k,dim], 0.95)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], 0.05)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], 0.5)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], 0.95)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(4,1:3)],qnew[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
}

for (dim in 1:3){
plotmcmc(qref[[dim]][,c(4,1:3)],qnew[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
}

x <- cbind(qref[[dim]]$iteration,qref[[dim]]$iteration,qref[[dim]]$iteration,qref[[dim]]$iteration,qref[[dim]]$iteration,qref[[dim]]$iteration)
y <- cbind(qref[[dim]]$V1,qref[[dim]]$V2,qref[[dim]]$V3,qnew[[dim]]$V1,qnew[[dim]]$V2,qnew[[dim]]$V3)

matplot(x,y,type="p")



#asymptotic variance
seed0 <- 39546
replicate <- 10
L_mcmc<-10000
final.ref <- 0
final.new <- 0
for (m in  1:replicate){
options_warfa<-list(seed=m*seed0,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta_ref
ref['individual'] <- m
final.ref <- rbind(final.ref,ref)

options_warfanew<-list(seed=m*seed0,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta
new['individual'] <- m
final.new <- rbind(final.new,new)
}





