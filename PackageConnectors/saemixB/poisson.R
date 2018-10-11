setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/PackageConnectors/saemixB/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 

  source('main.R')
  source('main_estep.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/PackageConnectors/saemixB")

library("mlxR")
library(MlxConnectors)
initializeMlxConnectors(software = "monolix")

################################################################ SAEMIX ####################################################################################################################################
project.file <- "mlxProjects/poisson/count1a_project.mlxtran"
loadProject(project.file)

getEstimatedPopulationParameters()
getEstimatedIndividualParameters()$saem
computePredictions(getEstimatedIndividualParameters()$saem)

getEstimatedIndividualParameters()$saem

count.data <- readDatamlx(project = project.file)
count.saemix <- count.data$Y

saemix.data<-saemixData(name.data=count.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("time","Y"),name.response=c("Y"), name.X="time")



countmodel<-function(psi,id,xidep) { 
  time<-xidep[,1]
  y<-xidep[,2]
  lambda<-psi[id,1]
  omega<-psi[id,2]
  logp <- -lambda + y*log(lambda) - factorial(y)
  return(logp)
}



saemix.model<-saemixModel(model=countmodel,description="count model",type="likelihood",   
  psi0=matrix(c(10,3),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),omega.init=matrix(c(0.3,0,0,0.3),ncol=2,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

##RUNS

K1 = 200
K2 = 50
iterations = 1:(K1+K2+1)
end = K1+K2

#Warfarin
options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,monolix=TRUE)
count<-saemix(saemix.model,saemix.data,options)

################################################################ SAEMIX ####################################################################################################################################


################################################################ MonolixProject ####################################################################################################################################
project.file <- "mlxProjects/warfarinmlx/warfarinPK_project.mlxtran"
loadProject(project.file)
warfa_data <- readDatamlx(project = project.file)
treat <- warfa_data$treatment[,c(1,3)]
warfarin.saemix <- merge(treat ,warfa_data$y_1,by="id")
warfarin.saemix <- warfarin.saemix[order(warfarin.saemix$id),]


saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y_1"), name.X="time")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  ypred <- 2
  return(ypred)
}

options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,nb.chains=1,monolix=TRUE)
warfa<-saemix(model=saemix.model,data=saemix.data,options_warfa)












################################################################ MonolixProject ####################################################################################################################################





