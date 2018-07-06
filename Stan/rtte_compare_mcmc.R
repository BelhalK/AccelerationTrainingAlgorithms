library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
library(rstan)
# save.image("rtte_mala.RData")
# save.image("rtte_mcmc_conv.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
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

load("rtte_mala.RData")
###RTTE
timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_new2/data/rtte_data.csv", header=T, sep=",")
# timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/data/rttellis.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
saemix.data_rtte<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
n <- length(unique(timetoevent.saemix$id))

timetoevent.model<-function(psi,id,xidep) {
T<-xidep[,1]
y<-xidep[,2]
N <- nrow(psi)
Nj <- length(T)
censoringtime = 20
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

saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))


##RUNS

# K1 = 200
# K2 = 100
# iterations = 1:(K1+K2+1)
# end = K1+K2

# #Weibull
# options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
# rtte <- cbind(iterations, rtte)


# options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5))
# rttenew<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenew))
# rtte <- cbind(iterations, rtte)



saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(10.17122,4.577724),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(0.3,0,0,0.3),ncol=2, 
  byrow=TRUE))


L_mcmc=10000
options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rtte)$eta_ref

options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rttenew)$eta


i <- 2
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F, av=0, sigma.val=0.01
  ,gamma.val=0.0001,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = i)
mala<-mcmc(saemix.model_rtte,saemix.data_rtte,options.mala)$eta



model <- 'data {
          int<lower=0> N_obs;// Number of observations
          vector[N_obs] time; //obs
          
          real<lower=0> lambda_pop;
          real<lower=0> beta_pop;
          real<lower=0> omega_lambda;
          real<lower=0> omega_beta;
        }
        parameters {
          vector<lower=0>[2] param;
        }

        model {
          //Priors
          param[1] ~ lognormal( lambda_pop , omega_lambda);
          param[2] ~ lognormal( beta_pop , omega_beta);
          time ~ weibull(param[1], param[2]);
        }'



modelstan <- stan_model(model_name = "rtte",model_code = model)

#NUTS using rstan
i <- 2
options.vi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(0,0,0,0,0,1),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
vi<-mcmc(saemix.model_rtte,saemix.data_rtte,options.vi)$eta




i <- 2
start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 2))


#quantiles
qlow <- 0.2
qmed <- 0.5
qhigh <- 0.8


qref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[[i]][1:k,dim], qlow)
    qref[[dim]][k,2] <- quantile(ref[[i]][1:k,dim], qmed)
    qref[[dim]][k,3] <- quantile(ref[[i]][1:k,dim], qhigh)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], qlow)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], qmed)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], qhigh)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(3,1:2)],qnew[[dim]][,c(3,1:2)],title=paste("quantiles",i,"dim", dim))
}

qmala <- list(mala[[i]][1:L_mcmc,],mala[[i]][1:L_mcmc,],mala[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qmala[[dim]][k,1] <- quantile(mala[[i]][1:k,dim], qlow)
    qmala[[dim]][k,2] <- quantile(mala[[i]][1:k,dim], qmed)
    qmala[[dim]][k,3] <- quantile(mala[[i]][1:k,dim], qhigh)
  }
  qmala[[dim]]$iteration <- 1:L_mcmc
}



iteration <- 1:L_mcmc
burn <- 100
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


q1mala <- data.frame(cbind(iteration,qmala[[1]][,1],qmala[[2]][,1]))
q2mala <- data.frame(cbind(iteration,qmala[[1]][,2],qmala[[2]][,2]))
q3mala <- data.frame(cbind(iteration,qmala[[1]][,3],qmala[[2]][,3]))
q1mala$quantile <- 1
q2mala$quantile <- 2
q3mala$quantile <- 3
quantmala <- rbind(q1mala[-c(1:burn),],q2mala[-c(1:burn),],q3mala[-c(1:burn),])



q1ref[,2] <- q1ref[,2] + 10 
q2ref[,2] <- q2ref[,2] + 10
q3ref[,2] <- q3ref[,2] + 10

q1new[,2] <- q1new[,2] + 10
q2new[,2] <- q2new[,2] + 10
q3new[,2] <- q3new[,2] + 10

q1mala[,2] <- q1mala[,2] + 10
q2mala[,2] <- q2mala[,2] + 10
q3mala[,2] <- q3mala[,2] + 10


q1ref[,3] <- q1ref[,3] + 3 
q2ref[,3] <- q2ref[,3] + 3
q3ref[,3] <- q3ref[,3] + 3

q1new[,3] <- q1new[,3] + 3
q2new[,3] <- q2new[,3] + 3
q3new[,3] <- q3new[,3] + 3

q1mala[,3] <- q1mala[,3] + 3
q2mala[,3] <- q2mala[,3] + 3
q3mala[,3] <- q3mala[,3] + 3


colnames(quantref) <- colnames(quantnew)<-colnames(quantmala)<-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")

plotquantile3(quantref,quantnew,quantmala)

