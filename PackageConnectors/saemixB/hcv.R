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

################################################################ MonolixProject ####################################################################################################################################

project.file <- "mlxProjects/hcv/hcv_project.mlxtran"
loadProject(project.file)

# getEstimatedPopulationParameters()
# computePredictions(getEstimatedIndividualParameters()$saem)
# computePredictions(getEstimatedIndividualParameters()$saem, individualIds = c(10,20))

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  s<-psi[id,1]
  d<-psi[id,2]
  beta<-psi[id,3]
  delta<-psi[id,4]
  p<-psi[id,5]
  c<-psi[id,6]
  eta<-psi[id,7]
  epsilon<-psi[id,8]
  ypred<-1
  return(ypred)
}

hcv_data <- readDatamlx(project = project.file)
treat <- hcv_data$treatment[,c(3)]
# hcv.saemix <- cbind(hcv_data$Y,treat)
hcv.saemix <- hcv_data$Y

saemix.data<-saemixData(name.data=hcv.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("time"),name.response=c("Y"), name.X="time")

saemix.model<-saemixModel(model=model1cpt,description="hcv",type="structural"
  ,psi0=matrix(c(1000,1,0.00005,0.05,20,5,0.9,0.7),ncol=8,byrow=TRUE, dimnames=list(NULL, c("s","d","beta","delta","p","c","eta","epsilon"))),fixed.estim=c(1,1,1,1,1,1,1,1),
  transform.par=c(1,1,1,1,1,1,1,1),omega.init=matrix(diag(8),ncol=8,byrow=TRUE),covariance.model=matrix(diag(8),ncol=8,byrow=TRUE))

K1 = 200
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,nb.chains=1,monolix=TRUE)
warfa<-saemix(saemix.model,saemix.data,options)



