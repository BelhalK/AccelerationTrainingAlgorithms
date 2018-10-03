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
# load("testttenuts.RData")
# save.image("testttenuts.RData")
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
beta <- psi[id,1]
lambda <- psi[id,2]
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

K1 = 200
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

# #Weibull
# options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
# rtte <- cbind(iterations, rtte)


# options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5))
# rttenew<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenew))
# rtte <- cbind(iterations, rtte)



saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(3,10),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(0,0),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

i <- 5
L_mcmc=1000

options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,
  L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1,
   nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, 
   map.range=c(0), indiv.index = i)
new<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rttenew)$eta

# model <- 'data {
#   int<lower=1> N_e; // Number of total observed events
#   int<lower=1> N_c; // Number of total censoring times 
#   vector<lower=0>[N_e] event_times; // Times of event occurrence
#   int<lower=0> cens_times; // Censoring times
#   real<lower=0> alpha_pop;
#   real<lower=0> sigma_pop;
#   real<lower=0> omega_alpha;
#   real<lower=0> omega_sigma;
# }

# parameters {
#   vector<lower=0>[2] beta;
# }

# model {
#   // prior
#   beta[1] ~ normal(alpha_pop, omega_alpha);
#   beta[2] ~ normal(sigma_pop, omega_sigma);
  
#   // likelihood
#   for (n_e in 1:N_e) {
#     target += weibull_lpdf(event_times[n_e] | beta[1], beta[2]) - 
#               weibull_lccdf(event_times[n_e] | beta[1], beta[2]);
#   } 

#   target += weibull_lccdf(cens_times | beta[1], beta[2]);

# }'


model <- 'data {
  int<lower=1> N_e; // Number of total observed events
  int<lower=1> N_c; // Number of total censoring time
  vector<lower=0>[N_e] event_times; // Times of event occurrence
  int<lower=0> cens_times; // Censoring time
  real<lower=0> beta_pop;
  real<lower=0> lambda_pop;
  real<lower=0> omega_beta;
  real<lower=0> omega_lambda;
}

parameters {
  vector<lower=0>[2] param;
}

model {
  // prior
  param[1] ~ normal(beta_pop, omega_beta);
  param[2] ~ normal(lambda_pop, omega_lambda);
  
  // likelihood
  target += weibull_lpdf(event_times | param[1], param[2]) - 
            weibull_lccdf(event_times | param[1], param[2]) +
            weibull_lccdf(cens_times | param[1], param[2]);
}'


modelstan <- stan_model(model_name = "rtte",model_code = model)

#NUTS using rstan
options.vi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
vi<-mcmc(saemix.model_rtte,saemix.data_rtte,options.vi)$eta


start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 2))

#quantiles
qlow <- 0.1
qmed <- 0.5
qhigh <- 0.9

qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], qlow)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], qmed)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], qhigh)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
}

qvi <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qvi[[dim]][k,1] <- quantile(vi[[i]][1:k,dim], qlow)
    qvi[[dim]][k,2] <- quantile(vi[[i]][1:k,dim], qmed)
    qvi[[dim]][k,3] <- quantile(vi[[i]][1:k,dim], qhigh)
  }
  qvi[[dim]]$iteration <- 1:L_mcmc
}



iteration <- 1:L_mcmc
burn <- 100

q1new <- data.frame(cbind(iteration,qnew[[1]][,1],qnew[[2]][,1]))
q2new <- data.frame(cbind(iteration,qnew[[1]][,2],qnew[[2]][,2]))
q3new <- data.frame(cbind(iteration,qnew[[1]][,3],qnew[[2]][,3]))
q1new$quantile <- 1
q2new$quantile <- 2
q3new$quantile <- 3
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])

q1vi <- data.frame(cbind(iteration,qvi[[1]][,1],qvi[[2]][,1]))
q2vi <- data.frame(cbind(iteration,qvi[[1]][,2],qvi[[2]][,2]))
q3vi <- data.frame(cbind(iteration,qvi[[1]][,3],qvi[[2]][,3]))
q1vi$quantile <- 1
q2vi$quantile <- 2
q3vi$quantile <- 3
quantnuts <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])

colnames(quantnew) <- colnames(quantnuts)<-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")


plotquantile3(quantnew,quantnew,quantnuts)
