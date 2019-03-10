

library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library("rCMA")
library(rstan)

  source('R/aaa_generics.R') 
  source('R/compute_LL.R') 
  source('R/func_aux.R') 
  source('R/func_distcond.R') 
  source('R/func_FIM.R')
  source('R/func_plots.R') 
  source('R/func_simulations.R') 
  source('R/main.R')
  source('R/main_estep.R')
  source('R/main_initialiseMainAlgo.R') 
  source('R/main_mstep.R') 
  source('R/SaemixData.R')
  source('R/SaemixModel.R') 
  source('R/SaemixRes.R') 
  # source('R/SaemixRes_c.R') 
  source('R/SaemixObject.R') 
  source('R/zzz.R') 
source('R/graphplot.R') 

bolus_data <- read.table("data/bolus1_data.txt", header=T)
saemix.data<-saemixData(name.data=bolus_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amt","time"),name.response=c("y"), name.X="time")


model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
  V<-psi[id,1]
  k<-psi[id,2]
  k12<-psi[id,3]
  k21<-psi[id,4]

  a0 <- k*k21*k31
  a1 <- k*k31 + k21*k31 + k21*k13 + k*k21 + k31*k12
  a2 <- k + k12 + k21
  p <- a1 - a2^2/3
  q <- 2*a2^3/27 - a1*a2/3 + a0
  r1 <- sqrt(-(p^3/27))
  r2 <- 2*r1^(1/3)
  phi <- acos(-q/(2*r1))/3
  alpha <- -(cos(phi)*r2 - a2/3)
  beta <- -(cos(phi + 2*pi/3)*r2 - a2/3)
  gamma <- -(cos(phi + 4*pi/3)*r2 - a2/3)

  A <- (k21 - alpha)*(k31 - alpha)/(V*(alpha-beta)*(alpha-gamma))
  B <- (k21 - beta)*(k31 - beta)/(V*(beta-alpha)*(beta-gamma))
  C <- (k21 - gamma)*(k31 - gamma)/(V*(gamma-beta)*(gamma-alpha))

  ypred<-dose*(A/alpha*(1-exp(-alpha*time))+ B/beta*(1-exp(-beta*time)) + C/gamma*(1-exp(-gamma*time)))
  return(ypred)
}
# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural",
  ,psi0=matrix(c(1,7),ncol=2,byrow=TRUE, dimnames=list(NULL, c("V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))


K1 = 100
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2


replicate = 3
seed0 = 395246
seed1 = 3952


final_rwm <- 0
final_mix <- 0
final_mala <- 0


# model <- 'data {
#           int<lower=0> N;// Number of observations
#           vector[N] time; //predictor
#           real dose; //predictor
#           vector[N] concentration;  //response
          
#           real beta1_pop;
#           real beta2_pop;
#           real beta3_pop;
#           real<lower=0> omega_beta1;
#           real<lower=0> omega_beta2;
#           real<lower=0> omega_beta3;
#           real<lower=0>  pres;
#         }
#         parameters {
#           vector<lower=0>[3] beta;
#         }
#         model {
#           //Priors
#           beta[1] ~ lognormal( beta1_pop , omega_beta1);
#           beta[2] ~ lognormal( beta2_pop , omega_beta2);
#           beta[3] ~ lognormal( beta3_pop , omega_beta3);

#           concentration ~ normal(dose*beta[1]/(beta[2]*(beta[1]-beta[3]))*(exp(-beta[3]*time)-exp(-beta[1]*time)), pres);
#         }'

# modelstan <- stan_model(model_name = "warfarin",model_code = model)
m=3
for (m in 1:replicate){
  print(m)
  l = list(c(50,2),c(80,3),c(60,4))
  # l = list(c(1,5,2,0,0,0),c(3,12,5,0,0,0),c(6,3,7,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(l[[m]],ncol=2,byrow=TRUE, dimnames=list(NULL, c("V","k"))),
  transform.par=c(1,1),omega.init=matrix(c(1/m,0,0,1/m),ncol=2,byrow=TRUE))

  options<-list(seed=seed0/m,map=F,fim=F,ll.is=T,nb.chains = 1,
   nbiter.mcmc = c(2,2,2,0,0,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0)
  bolus_ref<-data.frame(saemix(saemix.model,saemix.data,options)$par)
  bolus_ref <- cbind(iterations, bolus_ref)
  bolus_ref[,4:6] <- sqrt(bolus_ref[,4:6])
  bolus_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,bolus_ref[-1,])

  options.new<-list(seed=m*seed1,map=F,fim=F,ll.is=T,nb.chains = 1,
   nbiter.mcmc = c(2,2,2,2,0,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),
   map.range=c(1:4),nbiter.burn =0)
  bolus_new_ref<-data.frame(saemix(saemix.model,saemix.data,options.new)$par)
  bolus_mix <- cbind(iterations, bolus_new_ref)
  bolus_mix[,4:6] <- sqrt(bolus_mix[,4:6])
  bolus_mix['individual'] <- m
  final_mix <- rbind(final_mix,bolus_mix[-1,])

  #  options.mala<-list(seed=seed0/m,map=F,fim=F,ll.is=T,nb.chains = 1,
  #   nbiter.mcmc = c(2,2,2,0,2,0),
  #   nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(1),nbiter.burn =0,sigma.val=0.002,gamma.val=0.1)
  # bolus_mala_ref<-data.frame(saemix(saemix.model,saemix.data,options.mala)$par)
  # bolus_mala <- cbind(iterations, bolus_mala_ref)
  # bolus_mala[,4:6] <- sqrt(bolus_mala[,4:6])
  # bolus_mala['individual'] <- m
  # final_mala <- rbind(final_mala,bolus_mala[-1,])

  #  options.nuts<-list(seed=seed0/m,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,0,6),
  #   nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(1:15),nbiter.burn =0,sigma.val=0.002,gamma.val=0.1, modelstan = modelstan)
  # bolus_nuts_ref<-data.frame(saemix(saemix.model,saemix.data,options.nuts)$par)
  # bolus_nuts <- cbind(iterations, bolus_nuts_ref)
  # bolus_nuts[,5:7] <- sqrt(bolus_nuts[,5:7])
  # bolus_nuts['individual'] <- m
  # final_nuts <- rbind(final_nuts,bolus_nuts[-1,])

}


convpop <- graphConvMC_diffpk1(final_rwm[,c(1,3,7)],final_mix[,c(1,3,7)])
convvar <- graphConvMC_diffpk1(final_rwm[,c(1,5,7)],final_mix[,c(1,5,7)])

convpop <- graphConvMC_diffpk1(final_rwm[,c(1,2,7)],final_mix[,c(1,2,7)])
convvar <- graphConvMC_diffpk1(final_rwm[,c(1,4,7)],final_mix[,c(1,4,7)])
