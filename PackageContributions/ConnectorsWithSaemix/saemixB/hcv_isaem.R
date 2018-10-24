# save.image("Rdata/hcv.RData")
# load("Rdata/hcv.RData")
# save.image("Rdata/hcv_morechains.RData")
# load("Rdata/hcv_morechains.RData")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/PackageContributions/ConnectorsWithSaemix/saemixB/incrementalR")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('main.R')
  source('main_estep.R')
  source('main_estep_incremental.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  source('plots.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/PackageContributions/ConnectorsWithSaemix/saemixB")
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

project.file <- "mlxProjects/hcv/hcv_project.mlxtran"
loadProject(project.file)

# getEstimatedPopulationParameters()
# getEstimatedLogLikelihood()
# runLogLikelihoodEstimation(linearization = FALSE, wait = TRUE)
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

# saemix.model<-saemixModel(model=model1cpt,description="hcv",type="structural"
#   ,psi0=matrix(c(1000,1,0.00005,0.05,20,5,0.9,0.7),ncol=8,byrow=TRUE, dimnames=list(NULL, c("s","d","beta","delta","p","c","eta","epsilon"))),fixed.estim=c(1,1,1,1,1,1,1,1),
#   transform.par=c(1,1,1,1,1,1,3,3),omega.init=matrix(diag(8),ncol=8,byrow=TRUE),covariance.model=matrix(diag(8),ncol=8,byrow=TRUE))

# saemix.model<-saemixModel(model=model1cpt,description="hcv",type="structural"
#   ,psi0=matrix(c(1000,1,0.00005,0.05,20,5,0.9,0.7),ncol=8,byrow=TRUE, dimnames=list(NULL, c("s","d","beta","delta","p","c","eta","epsilon"))),fixed.estim=c(1,1,0,0,0,0,1,1),
#   transform.par=c(1,1,1,1,1,1,2,2),omega.init=matrix(diag(8),ncol=8,byrow=TRUE),covariance.model=matrix(diag(8),ncol=8,byrow=TRUE))


cov.model = matrix(0,nrow=8,ncol=8,byrow=TRUE)
cov.model[1,1] <- 1
cov.model[2,2] <- 1

saemix.model<-saemixModel(model=model1cpt,description="hcv",type="structural"
  ,psi0=matrix(c(1000,1,0.00005,0.05,20,5,0.9,0.7),ncol=8,byrow=TRUE,
   dimnames=list(NULL, c("s","d","beta","delta","p","c","eta","epsilon"))),
  fixed.estim=c(1,1,0,0,0,0,0,0),
  transform.par=c(1,1,1,1,1,1,2,2),omega.init=matrix(diag(8),ncol=8,byrow=TRUE),
  covariance.model=cov.model)



K1 = 2000
K2 = 1000
iterations = 1:(K1+K2)
end = K1+K2

runtime = 200
nchains = 50


options<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=100,sampling='randompass', duration = runtime)
hcv<-data.frame(saemix(saemix.model,saemix.data,options))
hcv <- cbind(iterations, hcv[-1,])
row_sub_ref  = apply(hcv, 1, function(row) all(row !=0 ))
hcv <- hcv[row_sub_ref,]
hcv$algo <- 'full'
hcv$iterations <- seq(0,runtime, length.out=length(hcv$iterations))



options25<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=25,sampling='randompass', duration = runtime)
hcv25<-data.frame(saemix(saemix.model,saemix.data,options25))
hcv25 <- cbind(iterations, hcv25[-1,])
row_sub_25  = apply(hcv25, 1, function(row) all(row !=0 ))
hcv25 <- hcv25[row_sub_25,]
hcv25$algo <- 'quarter'
hcv25$iterations <- seq(0,runtime, length.out=length(hcv25$iterations))




options50<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=50,sampling='randompass', duration = runtime)
hcv50<-data.frame(saemix(saemix.model,saemix.data,options50))
hcv50 <- cbind(iterations, hcv50[-1,])
row_sub_50  = apply(hcv50, 1, function(row) all(row !=0 ))
hcv50 <- hcv50[row_sub_50,]
hcv50$algo <- 'half'
hcv50$iterations <- seq(0,runtime, length.out=length(hcv50$iterations))



options75<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=75,sampling='randompass', duration = runtime)
hcv75<-data.frame(saemix(saemix.model,saemix.data,options75))
hcv75 <- cbind(iterations, hcv75[-1,])
row_sub_75  = apply(hcv75, 1, function(row) all(row !=0 ))
hcv75 <- hcv75[row_sub_75,]
hcv75$algo <- '75'
hcv75$iterations <- seq(0,runtime, length.out=length(hcv75$iterations))


options95<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=85,sampling='randompass', duration = runtime)
hcv95<-data.frame(saemix(saemix.model,saemix.data,options95))
hcv95 <- cbind(iterations, hcv95[-1,])
row_sub_95  = apply(hcv95, 1, function(row) all(row !=0 ))
hcv95 <- hcv95[row_sub_95,]
hcv95$algo <- '85'
hcv95$iterations <- seq(0,runtime, length.out=length(hcv95$iterations))







comparison <- 0
# comparison <- rbind(hcv,hcv)
# comparison <- rbind(hcv,hcv50)
comparison <- rbind(hcv,hcv25,hcv50)
# comparison <- rbind(hcv,hcv25,hcv50,hcv75,hcv95)
var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
prec <- seplot(var, title="comparison",legend=TRUE)


