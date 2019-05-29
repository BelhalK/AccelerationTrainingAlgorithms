setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incrementalwithMlxConnectors/R")
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
  # source('SaemixRes_c.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  
  source('main_incremental.R')
  source('main_estep_incremental.R')


  source('/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incrementalwithMlxConnectors/R/mixtureFunctions.R')
  source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incrementalwithMlxConnectors/plots.R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incrementalwithMlxConnectors")

library("rlist")
# library("mlxR")
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

library("mlxR")

project.file <- "hiv/hiv_project.mlxtran"
loadProject(project.file)

# getEstimatedPopulationParameters()
# computePredictions(getEstimatedIndividualParameters()$saem)
# computePredictions(getEstimatedIndividualParameters()$saem, individualIds = c(10,20))

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
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


hiv_data <- readDatamlx(project = project.file)
treat <- hiv_data$treatment[,c(3)]
# hiv.saemix <- cbind(hiv_data$Y,treat)
hiv.saemix <- hiv_data$Y


# warfa_data <- read.table("hiv/data/hiv1_data.txt", header=T)

saemix.data<-saemixData(name.data=hiv.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("time"),name.response=c("Y"), name.X="time")

# Default model, no covariate

saemix.model<-saemixModel(model=model1cpt,description="hiv",type="structural"
  ,psi0=matrix(c(1000,1,0.00005,0.05,20,5,0.9,0.7),ncol=8,byrow=TRUE, dimnames=list(NULL, c("s","d","beta","delta","p","c","eta","epsilon"))),fixed.estim=c(1,1,1,1,1,1,1,1),
  transform.par=c(1,1,1,1,1,1,1,1),omega.init=matrix(diag(8),ncol=8,byrow=TRUE),covariance.model=matrix(diag(8),ncol=8,byrow=TRUE))

# K1 = 2000
# K2 = 500
K1 = 300
K2 = 200
iterations = 1:(K1+K2)
end = K1+K2
batchsize25 = 25
batchsize50 = 50

seed0=39546


runtime <- 40
#####RWM#####
options<-list(seed=seed0,map=F,fim=F,ll.is=T,save.graphs=FALSE,nb.chains = 1,nbiter.mcmc = c(2,2,2,0),
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, 
 map.range=c(0), nb.replacement=100,sampling='randompass', duration = runtime)
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])
row_sub_ref  = apply(theo_ref, 1, function(row) all(row !=0 ))
theo_ref <- theo_ref[row_sub_ref,]
theo_ref$algo <- 'full'
theo_ref$iterations <- seq(0,10, length.out=length(theo_ref$iterations))



options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=T,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass', duration = runtime)
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])
row_sub25  = apply(theo_mix25, 1, function(row) all(row !=0 ))
theo_mix25 <- theo_mix25[row_sub25,]
theo_mix25$algo <- 'quarter'
theo_mix25$iterations <- seq(0,10, length.out=length(theo_mix25$iterations))



iternewkernel <- 1:2
runtime <- 60
#####NEW KERNEL#####
options<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1,
  nbiter.mcmc = c(2,2,2,2),
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, 
 map.range=c(iternewkernel), nb.replacement=100,sampling='randompass', duration = runtime)
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])
row_sub_ref  = apply(theo_ref, 1, function(row) all(row !=0 ))
theo_ref <- theo_ref[row_sub_ref,]
theo_ref$algo <- 'full'
theo_ref$iterations <- seq(0,runtime, length.out=length(theo_ref$iterations))



options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,2), nbiter.saemix = c(K1,K2),displayProgress=FALSE, 
  map.range=c(iternewkernel),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=50,sampling='randompass', duration = runtime)
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])
row_sub  = apply(theo_mix25, 1, function(row) all(row !=0 ))
theo_mix25 <- theo_mix25[row_sub,]
theo_mix25$algo <- 'quarter'
theo_mix25$iterations <- seq(0,runtime, length.out=length(theo_mix25$iterations))




comparison <- 0
comparison <- rbind(theo_ref[,],theo_mix25[,])
var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
prec <- seplot(var, title="comparison",legend=TRUE)
# assign(paste("prec", i, sep = ""), prec) 


