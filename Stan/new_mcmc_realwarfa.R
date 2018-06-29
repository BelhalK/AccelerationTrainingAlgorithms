library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
library(data.table)
library(rstan)
load("hmc_quantile.RData")
# save.image("hmc_quantile.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Stan/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('estep_mcmc.R')
  source('indiv_VI.R')
  source('variationalinferencelinear.R')
  source('main.R')
  source('main_estep.R')
  source('mcmc_final.R')
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

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Stan")



warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/Stan/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")


n <- length(unique(warfa_data$id))
model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]

  ypred<-dose*ka/(V*(ka-k))*(exp(-k*time)-exp(-ka*time))
  return(ypred)
}

# saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))


##RUNS

K1 = 400
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

# #Warfarin
# options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
# warfa <- cbind(iterations, warfa)


# options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(K1))
# warfanew<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew))
# warfanew <- cbind(iterations, warfanew)


# graphConvMC_twokernels(warfa,warfanew)


#compareMCMC

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.7,7.51,0.0178),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(0.2,0,0,0,0.18,0,0,0,0.03),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


L_mcmc=1000
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta_ref

options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta

#MALA
i <- 10
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F, av=0, sigma.val=0.01
  ,gamma.val=0.0001,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = i)
mala<-mcmc(saemix.model_warfa,saemix.data_warfa,options.mala)$eta


#RSTAN VB

model <- 'data {
          int<lower=0> N;// Number of observations
          vector[N] time; //predictor
          real dose; //predictor
          vector[N] concentration;  //response
          
          real beta1_pop;
          real beta2_pop;
          real beta3_pop;
          real<lower=0> omega_beta1;
          real<lower=0> omega_beta2;
          real<lower=0> omega_beta3;
          real<lower=0>  pres;
        }
        parameters {
          vector<lower=0>[3] beta;
        }
        model {
          //Priors
          beta[1] ~ lognormal( beta1_pop , omega_beta1);
          beta[2] ~ lognormal( beta2_pop , omega_beta2);
          beta[3] ~ lognormal( beta3_pop , omega_beta3);
          concentration ~ normal(dose*beta[1]/(beta[2]*(beta[1]-beta[3]))*(exp(-beta[3]*time)-exp(-beta[1]*time)), pres);
        }'



modelstan <- stan_model(model_name = "warfarin",model_code = model)

#NUTS using rstan
i <- 10
options.vi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(0,0,0,0,0,1),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
vi<-mcmc(saemix.model_warfa,saemix.data_warfa,options.vi)$eta

# #using the samples from vb (drawn from candidate KL posterior)
# options_warfavb<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,1,0,1),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), modelstan = modelstan)
# vb<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfavb)$eta

#Autocorrelation
rwm.obj <- as.mcmc(ref[[i]])
autocorr.plot(rwm.obj[,1]) + title("RWM SAEM Autocorrelation")

new.obj <- as.mcmc(new[[i]])
autocorr.plot(new.obj[,1]) + title("Laplace SAEM Autocorrelation")

# vb.obj <- as.mcmc(vb[[i]])
# autocorr.plot(vb.obj[,1]) + title("vb SAEM Autocorrelation")

vi.obj <- as.mcmc(vi[[i]])
autocorr.plot(vi.obj[,1]) + title("vb SAEM Autocorrelation")

#MSJD
mssd(ref[[i]][,1])
mssd(new[[i]][,1])
# mssd(mala[[i]][,1])
mssd(vi[[i]][,1])


start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 3))


#quantiles
qlow <- 0.2
qmed <- 0.5
qhigh <- 0.8


qref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[[i]][1:k,dim], qlow)
    qref[[dim]][k,2] <- quantile(ref[[i]][1:k,dim], qmed)
    qref[[dim]][k,3] <- quantile(ref[[i]][1:k,dim], qhigh)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], qlow)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], qmed)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], qhigh)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
}




# qvb <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
# for (dim in 1:3){
#   print(dim)
#   for (k in 1:L_mcmc){
#     qvb[[dim]][k,1] <- quantile(vb[[i]][1:k,dim], qlow)
#     qvb[[dim]][k,2] <- quantile(vb[[i]][1:k,dim], qmed)
#     qvb[[dim]][k,3] <- quantile(vb[[i]][1:k,dim], qhigh)
#   }
#   qvb[[dim]]$iteration <- 1:L_mcmc
# }



qvi <- list(vi[[i]][1:L_mcmc,],vi[[i]][1:L_mcmc,],vi[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qvi[[dim]][k,1] <- quantile(vi[[i]][1:k,dim], qlow)
    qvi[[dim]][k,2] <- quantile(vi[[i]][1:k,dim], qmed)
    qvi[[dim]][k,3] <- quantile(vi[[i]][1:k,dim], qhigh)
  }
  qvi[[dim]]$iteration <- 1:L_mcmc
}


qmala <- list(mala[[i]][1:L_mcmc,],mala[[i]][1:L_mcmc,],mala[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qmala[[dim]][k,1] <- quantile(mala[[i]][1:k,dim], qlow)
    qmala[[dim]][k,2] <- quantile(mala[[i]][1:k,dim], qmed)
    qmala[[dim]][k,3] <- quantile(mala[[i]][1:k,dim], qhigh)
  }
  qmala[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(4,1:3)],qnew2[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
}



iteration <- 1:L_mcmc
burn <- 100
q1ref <- data.frame(cbind(iteration,qref[[1]][,1],qref[[2]][,1],qref[[3]][,1]))
q2ref <- data.frame(cbind(iteration,qref[[1]][,2],qref[[2]][,2],qref[[3]][,2]))
q3ref <- data.frame(cbind(iteration,qref[[1]][,3],qref[[2]][,3],qref[[3]][,3]))
q1ref$quantile <- 1
q2ref$quantile <- 2
q3ref$quantile <- 3
quantref <- rbind(q1ref[-c(1:burn),],q2ref[-c(1:burn),],q3ref[-c(1:burn),])


q1new <- data.frame(cbind(iteration,qnew[[1]][,1],qnew[[2]][,1],qnew[[3]][,1]))
q2new <- data.frame(cbind(iteration,qnew[[1]][,2],qnew[[2]][,2],qnew[[3]][,2]))
q3new <- data.frame(cbind(iteration,qnew[[1]][,3],qnew[[2]][,3],qnew[[3]][,3]))
q1new$quantile <- 1
q2new$quantile <- 2
q3new$quantile <- 3
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])

colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
# plotquantile(quantref,quantnew)



# q1vb <- data.frame(cbind(iteration,qvb[[1]][,1],qvb[[2]][,1],qvb[[3]][,1]))
# q2vb <- data.frame(cbind(iteration,qvb[[1]][,2],qvb[[2]][,2],qvb[[3]][,2]))
# q3vb <- data.frame(cbind(iteration,qvb[[1]][,3],qvb[[2]][,3],qvb[[3]][,3]))
# q1vb$quantile <- 1
# q2vb$quantile <- 2
# q3vb$quantile <- 3
# quantvb <- rbind(q1vb[-c(1:burn),],q2vb[-c(1:burn),],q3vb[-c(1:burn),])


# colnames(quantvb)<-c("iteration","ka","V","k","quantile")


q1vi <- data.frame(cbind(iteration,qvi[[1]][,1],qvi[[2]][,1],qvi[[3]][,1]))
q2vi <- data.frame(cbind(iteration,qvi[[1]][,2],qvi[[2]][,2],qvi[[3]][,2]))
q3vi <- data.frame(cbind(iteration,qvi[[1]][,3],qvi[[2]][,3],qvi[[3]][,3]))
q1vi$quantile <- 1
q2vi$quantile <- 2
q3vi$quantile <- 3
quantvi <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])
colnames(quantvi)<-c("iteration","ka","V","k","quantile")



q1mala <- data.frame(cbind(iteration,qmala[[1]][,1],qmala[[2]][,1],qmala[[3]][,1]))
q2mala <- data.frame(cbind(iteration,qmala[[1]][,2],qmala[[2]][,2],qmala[[3]][,2]))
q3mala <- data.frame(cbind(iteration,qmala[[1]][,3],qmala[[2]][,3],qmala[[3]][,3]))
q1mala$quantile <- 1
q2mala$quantile <- 2
q3mala$quantile <- 3
quantmala <- rbind(q1mala[-c(1:burn),],q2mala[-c(1:burn),],q3mala[-c(1:burn),])


colnames(quantmala)<-c("iteration","ka","V","k","quantile")


plotquantile3(quantnew,quantnew,quantvi)
# plotquantile3(quantref,quantnew,quantvb)
plotquantile3(quantref,quantnew,quantvi)
plotquantile3(quantref,quantnew,quantmala)

plotquantile4(quantref,quantnew,quantvi,quantmala)


# gelman.plot(mcmc.list(as.mcmc(ref[[10]])), bin.width = 10, max.bins = 50,confidence = qhigh, transform = FALSE, autoburnin=TRUE, auto.layout = TRUE)
# geweke.plot(mcmc.list(as.mcmc(ref[[10]])), frac1=0.1, frac2=0.5)
# geweke.plot(mcmc.list(as.mcmc(new[[10]])), frac1=0.1, frac2=0.5)
