## Output states of the MCMC for all individuals of the population
library(saemix)
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
options.mcmc<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0),nb.chains=1)
states.ref<-mcmc(saemix.model,saemix.data,options.mcmc)$eta

options.mcmc.new<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6),nb.chains=1)
states.new<-mcmc(saemix.model,saemix.data,options.mcmc.new)$eta



#### with the MLE #####
  ###SAEM###
K1 = 400
K2 = 100

options.saem<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0)
parpop<-saemix(saemix.model,saemix.data,options.saem)

omega.estimate <- parpop["results"]["omega"]
res.estimate <- parpop["results"]["respar"]
parpop.estimate <- parpop["results"]["fixed.effects"]


  ###MCMC###
saemix.model.mcmc<-saemixModel(model=model1cpt,description="warfarin",type="structural",psi0=matrix(parpop.estimate,ncol=length(parpop.estimate),byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),type="structural",
  transform.par=c(1,1,1),omega.init=omega.estimate,
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

L_mcmc=100
options.mcmc<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0),nb.chains=1)
sates<-mcmc(saemix.model.mcmc,saemix.data,options.mcmc)$eta
####################

