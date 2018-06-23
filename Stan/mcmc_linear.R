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
# save.image("vb_linear.RData")
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

model <- inlineModel("


[INDIVIDUAL]
input = {base_pop, slope_pop,  omega_base, omega_slope}
DEFINITION:
base = {distribution=normal, reference=base_pop, sd=omega_base}
slope  = {distribution=normal, reference=slope_pop,  sd=omega_slope }

[LONGITUDINAL]
input = {base, slope,a}
EQUATION:
C = base + slope*t
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=100

param   <- c(
  base_pop  = 140,    omega_base  = 0.5,
  slope_pop   = 1,   omega_slope   = 0.4, a =1)
  
res <- simulx(model     = model,
              parameter = param,
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(1,10,by=1)))

data<- res$y
saemix.data<-saemixData(name.data=data,header=TRUE,
  name.group=c("id"),name.predictors=c("time"),name.response=c("y"),
  units=list(x="yr",y="cm"))

saemix.model<-saemixModel(model=growth.linear,description="Linear model",type="structural",
  psi0=matrix(c(130,2),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2,byrow=TRUE),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE), 
  error.model="constant")


# #Warfarin
# options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# linear<-data.frame(saemix(saemix.model,saemix.data,options))
# linear <- cbind(iterations, linear)



saemix.model<-saemixModel(model=growth.linear,description="Linear model",type="structural",
  psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(0,0),covariance.model=matrix(c(0.2,0,0,0.2),ncol=2,byrow=TRUE),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE), 
  error.model="constant")
L_mcmc=1000
#RWM mcmc
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model,saemix.data,options_warfa)$eta_ref

#New kernel mcmc
options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model,saemix.data,options_warfanew)$eta


#RSTAN 

model <- 'data {
          int<lower=0> N;// Number of observations
          vector[N] age; //predictor
          vector[N] height;  //response
          
          real beta1_pop;
          real beta2_pop;
          real<lower=0> omega_beta1;
          real<lower=0> omega_beta2;
          real<lower=0>  pres;
        }
        parameters {
          vector[2] beta;
        }
        model {
          //Priors
          beta[1] ~ normal( beta1_pop , omega_beta1);
          beta[2] ~ normal( beta2_pop , omega_beta2);
          height ~ normal(beta[1] + beta[2] * age, 1);
        }'


modelstan <- stan_model(model_name = "oxboys",model_code = model)

#NUTS using rstan
options.vi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,0,1),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), modelstan = modelstan)
vi<-mcmc(saemix.model,saemix.data,options.vi)$eta



i <- 2
#Autocorrelation
rwm.obj <- as.mcmc(ref[[i]])
autocorr.plot(rwm.obj[,1]) + title("RWM SAEM Autocorrelation")

new.obj <- as.mcmc(new[[i]])
autocorr.plot(new.obj[,1]) + title("Laplace SAEM Autocorrelation")


vi.obj <- as.mcmc(vi[[i]])
autocorr.plot(vi.obj[,1]) + title("VI SAEM Autocorrelation")

#MSJD
mssd(ref[[i]][,1])
mssd(new[[i]][,1])
mssd(vi[[i]][,1])


#quantiles

i <- 2
qref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[[i]][1:k,dim], 0.05)
    qref[[dim]][k,2] <- quantile(ref[[i]][1:k,dim], 0.5)
    qref[[dim]][k,3] <- quantile(ref[[i]][1:k,dim], 0.95)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], 0.05)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], 0.5)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], 0.95)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(4,1:3)],qnew[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
}


qvi <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qvi[[dim]][k,1] <- quantile(vi[[i]][1:k,dim], 0.05)
    qvi[[dim]][k,2] <- quantile(vi[[i]][1:k,dim], 0.5)
    qvi[[dim]][k,3] <- quantile(vi[[i]][1:k,dim], 0.95)
  }
  qvi[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(4,1:3)],qnew2[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
}



iteration <- 1:L_mcmc
burn <- 1
q1ref <- data.frame(cbind(iteration,qref[[1]][,1],qref[[2]][,1]))
q2ref <- data.frame(cbind(iteration,qref[[1]][,2],qref[[2]][,2]))
q3ref <- data.frame(cbind(iteration,qref[[1]][,3],qref[[2]][,3]))
q1ref$quantile <- 1
q2ref$quantile <- 2
q3ref$quantile <- 3
quantref <- rbind(q1ref[-c(1:burn),],q2ref[-c(1:burn),],q3ref[-c(1:burn),])


q1new <- data.frame(cbind(iteration,qnew[[1]][,1],qnew[[2]][,1]))
q2new <- data.frame(cbind(iteration,qnew[[1]][,2],qnew[[2]][,2]))
q3new <- data.frame(cbind(iteration,qnew[[1]][,3],qnew[[2]][,3]))
q1new$quantile <- 1
q2new$quantile <- 2
q3new$quantile <- 3
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])


colnames(quantref) <- colnames(quantnew)<-c("iteration","A","B","quantile")



q1vi <- data.frame(cbind(iteration,qvi[[1]][,1],qvi[[2]][,1]))
q2vi <- data.frame(cbind(iteration,qvi[[1]][,2],qvi[[2]][,2]))
q3vi <- data.frame(cbind(iteration,qvi[[1]][,3],qvi[[2]][,3]))
q1vi$quantile <- 1
q2vi$quantile <- 2
q3vi$quantile <- 3
quantvi <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])


colnames(quantvi)<-c("iteration","A","B","quantile")


plotquantile3(quantref,quantnew,quantvi)

