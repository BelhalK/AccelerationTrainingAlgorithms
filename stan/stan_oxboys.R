
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/cmaes/Dir")
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_cov.R') 
  source('func_distcond.R') 
  source('func_FIM.R') 
  source('func_ggplot2.R') 
  source('func_plots.R') 
  source('func_simulations.R') 
  source('ggplot2_global.R') 
  # source('KL.R') 
  #source('vi.R') 
  source('global.R')
  source('main.R')
  source('mcmc_main.R') 
  # source('main_estep.R')
  source('main_estep_mcmc.R') 
  source('main_estep_morekernels.R') 
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('plots_ggplot2.R') 
  source('saemix-package.R') 
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  source('main_new.R')
  source('main_estep_new2.R')
  source('main_new_mix.R')
  source('main_estep_mix.R')
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/stan/")
source("mixtureFunctions.R")
source('main_estep_stan.R')
source('main_stan.R')

library(lme4)
library(rstan)
library(shinystan)#for great model viz
library(ggplot2)#for great viz in general
data(Penicillin)

library(RColorBrewer)
library("mlxR")
library("psych")
library("coda")
library("Matrix")
#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

# Doc
oxboys.saemix<-read.table( "oxboys.saemix.tab",header=T,na=".")
saemix.data<-saemixData(name.data=oxboys.saemix,header=TRUE,
  name.group=c("Subject"),name.predictors=c("age"),name.response=c("height"),
  units=list(x="yr",y="cm"))


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



saemix.model<-saemixModel(model=growth.linear,description="Linear model",
  psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(1,0),covariance.model=matrix(c(1,1,1,1),ncol=2,byrow=TRUE), 
  error.model="constant")


K1 = 10
K2 = 2
iterations = 1:(K1+K2+1)
gd_step = 0.01
end = K1+K2



options.ref<-list(seed=395246,map=F,fim=F,ll.is=F,displayProgress=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0),nbiter.saemix = c(K1,K2),nbiter.burn =0)
ox_ref<-data.frame(saemix_stan(saemix.model,saemix.data,options.ref))
ox_ref <- cbind(iterations, ox_ref)
ox_ref[end,]
graphConvMC_twokernels(ox_ref,ox_ref, title="new kernel")

modelstan <- 'data {
          // First we declare all of our variables in the data block
          int<lower=0> N;// Number of observations
          vector[N] age; //Identify our predictor as a vector
          vector[N] height;  //Identify our outcome variable as a vector
        }
        parameters {
          vector[2] beta; //Our betas are a vector of length 2 (intercept and slope)
          real<lower=0> sigma; //error parameter
        }
        model {
          //Priors
          beta[1] ~ normal( 1 , .001); //intercept
          beta[2] ~ normal( 0 , 0.1 ); //slope
          sigma ~ uniform( 0 , 1 ); //error
          height ~ normal(beta[1] + beta[2] * age, sigma);
        }'

options.stan<-list(seed=395246,map=F,fim=F,ll.is=F,displayProgress=FALSE,nb.chains = 1, nbiter.mcmc = c(0,0,0,1),nbiter.saemix = c(K1,K2),nbiter.burn =0, modelstan = modelstan)
ox_stan<-data.frame(saemix_stan(saemix.model,saemix.data,options.stan))
ox_stan <- cbind(iterations, ox_stan)
graphConvMC_twokernels(ox_stan,ox_stan, title="new kernel")

graphConvMC_twokernels(ox_ref,ox_stan, title="new kernel")

