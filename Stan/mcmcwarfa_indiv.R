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
load("RData/hmc_quantile_indiv.RData")
# save.image("hmc_quantile_indiv.RData")
# save.image("hmc_quantile_indiv_student.RData")
save.image("hmc_quantile_indiv_averagechains.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
source('R/aaa_generics.R') 
source('R/compute_LL.R') 
source('R/func_aux.R') 
source('R/func_distcond.R') 
source('R/func_FIM.R')
source('R/func_plots.R') 
source('R/func_simulations.R') 
source('R/estep_mcmc.R')
source('R/indiv_VI.R')
source('R/variationalinferencelinear.R')
source('R/main.R')
source('R/main_estep.R')
source('R/mcmc_final.R')
source('R/main_initialiseMainAlgo.R') 
source('R/main_mstep.R') 
source('R/check_linearvslaplace.R')
source('R/SaemixData.R')
source('R/SaemixModel.R') 
source('R/SaemixRes.R') 
# source('R/SaemixRes_c.R') 
source('R/SaemixObject.R') 
source('R/zzz.R')
source('R/graphplot.R')





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

# saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(0.7,7.51,0.0178),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(0.5,0,0,0,0.2,0,0,0,0.03),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,8,0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(0.2,0,0,0,0.18,0,0,0,0.03),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


L_mcmc=20000
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=i)
# ref<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=i)
# new<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta

new.student<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta

#MALA

i <- 10
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F, av=0, sigma.val=0.002
  ,gamma.val=0.1,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = i)
mala<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.mala)$eta


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
L_mcmc <- 100000
i <- 10
options.vi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
vi<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.vi)$eta


#ADVI for VI post outputs
#Calculate mu and gamma of ELBO optimization
variational.post.options<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nb.chains=1,
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),
  modelstan = modelstan, indiv.index = i)

variational.post<-indiv.variational.inference(saemix.model_warfa,saemix.data_warfa,variational.post.options)
mu.vi <- variational.post$mu
Gamma.vi <- variational.post$Gamma
etamap <- variational.post$map
Gammamap <- variational.post$Gammamap

# #using the output of ADVI (drawn from candidate KL posterior)
# test <- etamap
# # test[i,] <- etamap[i,] +0.01
# test[i,] <- mu.vi
eta.vi <- etamap
Gammavi <- Gammamap
# eta.vi[i,] <- mu.vi
# Gammavi[[i]] <- Gamma.vi
options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=eta.vi,Gamma = Gammavi,
        nbiter.mcmc = c(0,0,0,0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),
        nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index = i)
advi<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfavi)$eta

abs(mu.vi - etamap[i,])/abs(etamap[i,])
norm(Gamma.vi - Gammamap[[i]])/norm(Gammamap[[i]])
abs(Gamma.vi - Gammamap[[i]])/abs(Gammamap[[i]])
(Gamma.vi - Gammamap[[i]])/(Gammamap[[i]])



#Autocorrelation
rwm.obj <- as.mcmc(ref)
refac <- autocorr.plot(rwm.obj[,1]) 

new.obj <- as.mcmc(new)
newac <- autocorr.plot(new.obj[,1]) 

mala.obj <- as.mcmc(mala)
malaac <- autocorr.plot(mala.obj[,1]) 

vi.obj <- as.mcmc(vi)
viac <- autocorr.plot(vi.obj[,1]) 


advi.obj <- as.mcmc(advi)
adviac <- autocorr.plot(advi.obj[,1]) 

#Autocorrelation
par(mfrow=c(1,5))
acf(ref[,1], main="RWM")
acf(new[,1], main="IMH (Gaussian)")
acf(new.student[,1], main="IMH (Student)")
acf(mala[,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(advi[,1], main="ADVI")

par(mfrow=c(1,5))
acf(ref[,2], main="RWM")
acf(new[,2], main="IMH")
acf(mala[,2], main="MALA")
acf(vi[,2], main="NUTS")
acf(advi[,2], main="ADVI")


par(mfrow=c(1,5))
acf(ref[,3], main="RWM")
acf(new[,3], main="IMH")
acf(mala[,3], main="MALA")
acf(vi[,3], main="NUTS")
acf(advi[,3], main="ADVI")


#MSJD
mssd(ref[,1])
mssd(new[,1])
mssd(new.student[,1])
mssd(mala[,1])
mssd(advi[,1])
mssd(vi[,1])



# #Autocorrelation
# rwm.obj <- as.mcmc(states.ref[[10]])
# autocorr.plot(rwm.obj[,1]) + title("RWM Autocorrelation")

# new.obj <- as.mcmc(states.new[[10]])
# autocorr.plot(new.obj[,1]) + title("Laplace Autocorrelation")

# mala.obj <- as.mcmc(states.mala[[10]])
# autocorr.plot(mala.obj[,1]) + title("MALA Autocorrelation")

# nuts.obj <- as.mcmc(states.nuts[[10]])
# autocorr.plot(nuts.obj[,1]) + title("NUTS Autocorrelation")





L_mcmc=200
i <- 10
nchains <- 5
#REF
listofrefchains <- list(ref,ref)
for (m in 1:nchains){
  options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=i)
  ref<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta
  listofrefchains[[m]] <- ref
}


#new
listofnewchains <- list(ref,ref)
for (m in 1:nchains){
options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=i)
new<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta
listofnewchains[[m]] <- new
}
# new.student<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta



#MALA
listofmalachains <- list(ref,ref)
for (m in 1:nchains){
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F, av=0, sigma.val=0.002
  ,gamma.val=0.1,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = i)
mala<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.mala)$eta
listofmalachains[[m]] <- mala
}


#NUTS
listofnutschains <- list(ref,ref)
# listofnutschains <- 0
for (m in 1:nchains){
options.vi<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
vi<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.vi)$eta
  # listofnutschains <- listofnutschains + nuts
listofnutschains[[m]] <- vi
}



#ADVI
listofadvichains <- list(ref,ref)
for (m in 1:nchains){
options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=eta.vi,Gamma = Gammavi,
        nbiter.mcmc = c(0,0,0,0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),
        nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index = i)
advi<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfavi)$eta
listofadvichains[[m]] <- advi
}

