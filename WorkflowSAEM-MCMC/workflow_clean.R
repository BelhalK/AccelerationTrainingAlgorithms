library(saemix)
source('R/mcmc.R')  
# save.image("outputsaem.RData")
# load("outputsaem.RData")

warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/WorkflowSAEM-MCMC/data/warfarin_data.txt", header=T)
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


saemix.model<-saemixModel(model=model1cpt,description="warfarin",psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

K1 = 400
K2 = 100

#SAEM estimation
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0)
warfa<-saemix(saemix.model,saemix.data,options_warfa)

# parpop <- warfa["results"]["parpop"]
omega.estimate <- warfa["results"]["omega"]
res.estimate <- warfa["results"]["respar"]
parpop.estimate <- warfa["results"]["fixed.effects"]


#compareMCMC
saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural",psi0=matrix(parpop.estimate,ncol=length(parpop.estimate),byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=omega.estimate,
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

L_mcmc=100
options_warfa<-list(seed=39546,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0))
ref<-mcmc(saemix.model,saemix.data,options_warfa)$eta
