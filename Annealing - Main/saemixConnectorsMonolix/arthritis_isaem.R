# save.image("arthritis_40chains.RData")
# load("arthritis_40chains.RData")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Annealing - Main/saemixConnectorsMonolix/R")
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
  source("plots.R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Annealing - Main/saemixConnectorsMonolix")
library("rlist")
library("mlxR")
library(MlxConnectors)
initializeMlxConnectors(software = "monolix")

library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
require(madness)
################################################################ MonolixProject ####################################################################################################################################

project.file <- "mlxProjects/arthritis/arthritis_projet.mlxtran"
loadProject(project.file)

getEstimatedPopulationParameters()
getEstimatedLogLikelihood()
runLogLikelihoodEstimation(linearization = FALSE, wait = TRUE)
computePredictions(getEstimatedIndividualParameters()$saem)
computePredictions(getEstimatedIndividualParameters()$saem, individualIds = c(10,20))

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tauD<-psi[id,1]
  tauS<-psi[id,2]
  ypred<-1
  return(ypred)
}

arthritis_data <- readDatamlx(project = project.file)
treat <- arthritis_data$treatment[,c(3)]
# arthritis.saemix <- cbind(arthritis_data$Y,treat)
arthritis.saemix <- arthritis_data$y

saemix.data<-saemixData(name.data=arthritis.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("time"),name.response=c("y"), name.X="time")

saemix.model<-saemixModel(model=model1cpt,description="arthritis",type="structural"
  ,psi0=matrix(c(10,10),ncol=2,byrow=TRUE, dimnames=list(NULL, c("tauD","tauS"))),
  fixed.estim=c(1,1),transform.par=c(1,1),omega.init=matrix(diag(2),ncol=2,byrow=TRUE),
  covariance.model=matrix(diag(2),ncol=2,byrow=TRUE), error.model = "proportional")

K1 = 1000
K2 = 500
iterations = 1:(K1+K2)
end = K1+K2

runtime = 10
nchains = 1

options<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=100,sampling='randompass', duration = runtime)

arthritis<-data.frame(saemix(saemix.model,saemix.data,options))
arthritis <- cbind(iterations, arthritis[-1,])
row_sub_ref  = apply(arthritis, 1, function(row) all(row !=0 ))
arthritis <- arthritis[row_sub_ref,]
arthritis$algo <- 'full'
arthritis$iterations <- seq(0,runtime, length.out=length(arthritis$iterations))





K1 = 100
K2 = 50
iterations = 1:(K1+K2+1)
end = K1+K2
nchains = 1

replicate = 3
seed0 = 395246

#RWM
final_optim <- 0
final_av <- 0
final_avnew <- 0
final_bayes <- 0

final_annealing <- 0


for (m in 1:replicate){
  print(m)
  l = list(c(8,8),c(14,14),c(10,10),c(8,8))
  

saemix.model_arthritisnovar<-saemixModel(model=model1cpt,description="arthritis",type="structural"
  ,psi0=matrix(l[[m]],ncol=2,byrow=TRUE, dimnames=list(NULL, c("tauD","tauS"))),
  fixed.estim=c(1,1),transform.par=c(1,1),omega.init=matrix(diag(2),ncol=2,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,1),ncol=2,byrow=TRUE), error.model = "proportional")


  #No var
  ##### Optim (fmin search)
  options_arthritis_with<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
   duration = runtime,an=FALSE,coeff=1,nbiter.burn =0, av=0)
  arthritis_optim<-data.frame(saemix(saemix.model_arthritisnovar,saemix.data,options_arthritis_with))
  arthritis_optim <- cbind(iterations, arthritis_optim)
  arthritis_optim['individual'] <- m
  final_optim <- rbind(final_optim,arthritis_optim)
  
  #### AV
  options_arthritis_with<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=K1,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
   duration = runtime,an=FALSE,coeff=1,nbiter.burn =0, av=0)
  arthritis_withav<-data.frame(saemix(saemix.model_arthritisnovar,saemix.data,options_arthritis_with))
  arthritis_withav <- cbind(iterations, arthritis_withav)
  arthritis_withav['individual'] <- m
  final_av <- rbind(final_av,arthritis_withav)

  # ##### AV and new kernel
  # options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6,0),
  # nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(1), av=1)
  # arthritis_newkernelav<-data.frame(saemix(saemix.model_arthritisnovar,saemix.data,options_newkernel))
  # arthritis_newkernelav <- cbind(iterations, arthritis_newkernelav)
  # arthritis_newkernelav['individual'] <- m
  # final_avnew <- rbind(final_avnew,arthritis_newkernelav)


  # ##### Annealing MCMC
  options_an<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
   duration = runtime, an=TRUE,coeff=0.005,
  map.range=c(0), av=1)
  arthritis_an<-data.frame(saemix(saemix.model_arthritisnovar,saemix.data,options_an))
  arthritis_an <- cbind(iterations, arthritis_an)
  arthritis_an['individual'] <- m
  final_annealing<- rbind(final_annealing,arthritis_an)

  ##### pseudo bayesian

  options_arthritis_with<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
   duration = runtime,an=FALSE,coeff=1, 
    nbiter.burn =0, av=0,map.range=c(0))
  arthritis_bayes<-data.frame(saemix(saemix.model_arthritisnovar,saemix.data,options_arthritis_with))
  arthritis_bayes <- cbind(iterations, arthritis_bayes)
  arthritis_bayes['individual'] <- m
  final_bayes <- rbind(final_bayes,arthritis_bayes)

}


graphConvMC_diff4(final_optim,final_av,final_annealing,final_bayes, title="")
#black: optim
#blue: av
#red: annealing
#green: bayes


