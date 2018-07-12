## Output states of the MCMC for all individuals of the population
library(saemix)
library(rstan)
source('R/mcmc.R')  
source('R/func_aux.R') 
source('R/main_initialiseMainAlgo.R') 
source('R/SaemixModel.R') 
source('R/SaemixObject.R') 

warfa_data <- read.table("data/warfarin_data.txt", header=T)
saemix.data<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")


model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]

  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}


saemix.model<-saemixModel(model=model1cpt,description="warfarin",psi0=matrix(c(1,7,1,0,0,0),
  ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),type="structural",
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

L_mcmc=100
options.mcmc<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0,0,0),nb.chains=1)
states.ref<-mcmc(saemix.model,saemix.data,options.mcmc)$eta

options.mcmc.new<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0),nb.chains=1)
states.new<-mcmc(saemix.model,saemix.data,options.mcmc.new)$eta


options.mcmc.mala<-list(seed=39546,L_mcmc=L_mcmc,gamma.val=0.1,nbiter.mcmc = c(0,0,0,0,6,0),nb.chains=1)
states.mala<-mcmc(saemix.model,saemix.data,options.mcmc.mala)$eta


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
i <- 1
options.mcmc.nuts<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,0,1),nb.chains=1,
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i))
states.nuts<-mcmc(saemix.model,saemix.data,options.mcmc.nuts)$eta

