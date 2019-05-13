## Output states of the MCMC for all individuals of the population
library(saemix)
library(rstan)
source('R/mcmc.R')  
source('R/func_aux.R') 
source('R/main_initialiseMainAlgo.R') 
source('R/SaemixModel.R') 
source('R/SaemixObject.R') 

timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_new2/data/rtte_data.csv", header=T, sep=",")
# timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/data/rttellis.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
saemix.data<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
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



saemix.model<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(10,3),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),omega.init=matrix(c(0.3,0,0,0.3),ncol=2,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))


L_mcmc=100
options.mcmc<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0,0,0),nb.chains=1)
states.ref<-mcmc(saemix.model,saemix.data,options.mcmc)$eta

options.mcmc.new<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0),nb.chains=1)
states.new<-mcmc(saemix.model,saemix.data,options.mcmc.new)$eta


options.mcmc.mala<-list(seed=39546,L_mcmc=L_mcmc,gamma.val=0.0001,nbiter.mcmc = c(0,0,0,0,6,0),nb.chains=1)
states.mala<-mcmc(saemix.model,saemix.data,options.mcmc.mala)$eta


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
  param[2] ~ lognormal(beta_pop, omega_beta);
  param[1] ~ lognormal(lambda_pop, omega_lambda);
  
  // likelihood
  target += weibull_lpdf(event_times | param[2], param[1]) - 
            weibull_lccdf(event_times | param[2], param[1]) +
            weibull_lccdf(cens_times | param[2], param[1]);
}'


modelstan <- stan_model(model_name = "rtte",model_code = model)
options.mcmc.nuts<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,0,1),nb.chains=1,
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i))
states.nuts<-mcmc(saemix.model,saemix.data,options.mcmc.nuts)$eta

