setwd("/Users/karimimohammedbelhal/Desktop/averageiterate")
source('mixtureFunctions.R') 

setwd("/Users/karimimohammedbelhal/Desktop/averageiterate/avg")
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
  source('main_initialiseMainAlgoavg.R')
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 

library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
###RTTE
timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/paramToRV/data/rtte_data.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
saemix.data_rtte<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
timetoevent.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2]
  N <- nrow(psi)
  Nj <- length(T)
  censoringtime = 20
  lambda <- psi[id,1]
  beta <- psi[id,2]
  init <- which(T==0)
  cens <- which(T==censoringtime)
  ind <- setdiff(1:Nj, append(init,cens))
  hazard <- (beta/lambda)*(T/lambda)^(beta-1)
  H <- (T/lambda)^beta
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens] + H[cens-1]
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
  return(logpdf)
}

saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

##RUNS

K1 = 500
K2 = 500
iterations = 1:(K1+K2+1)
end = K1+K2


#With var no sa
options.ref<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=0)
rtte.ref<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options.ref))
rtte.ref <- cbind(iterations, rtte.ref)
# graphConvMC_twokernels(rtte.ref,rtte.ref)

options.avg<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=1)
rtte.avg<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options.avg))
rtte.avg <- cbind(iterations, rtte.avg)

# graphConvMC_twokernels(rtte.avg,rtte.avg)
graphConvMC_twokernels(rtte.ref,rtte.avg)

graphConvMC_twokernels(rtte.ref[K1:end,],rtte.avg[K1:end,])

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(1:5), av=0,avg=0)
rtte_newkernel<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_newkernel))
rtte_newkernel <- cbind(iterations, rtte_newkernel)

options_newavg<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(1:5), av=0,avg=1)
rtte_newavg<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_newavg))
rtte_newavg <- cbind(iterations, rtte_newavg)
graphConvMC_twokernels(rtte_without,rtte_newkernel)



final_ref <- rtte.ref
final_avg <- rtte.avg

seed0 = 39546
replicate = 10

for (j in 1:replicate){

#   model2 <- inlineModel("

#   [LONGITUDINAL]
#   input = {beta,lambda}  

#   EQUATION:
#   h=(beta/lambda)*(t/lambda)^(beta-1)

#   DEFINITION:
#   e = {type               = event, 
#        rightCensoringTime = 6,  
#        hazard             = h}
#   [INDIVIDUAL]
#   input={lambda_pop, o_lambda,beta_pop, o_beta}
                        
#   DEFINITION:
#   lambda  ={distribution=lognormal, prediction=lambda_pop,  sd=o_lambda}
#   beta  ={distribution=lognormal, prediction=beta_pop,  sd=o_beta}
#        ")


#   p <- c(lambda_pop=10, o_lambda=1,
#          beta_pop = 2,o_beta = 1)
#   h <- list(name='h', time=seq(0, 6, by=1))
#   e <- list(name='e', time=0)

#   N <- 10
#   res <- simulx(model     = model2, 
#                 settings  = list(seed=123),
#                 parameter = p, 
#                 output    = list(h, e), 
#                  group     = list(size = N))



# writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Desktop/paramToRV/data/rtte2.csv")
# timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/paramToRV/data/rtte2.csv", header=T, sep=",")
# timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
  saemix.data<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
  saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))
  print(j)

  options.ref<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=0)
  theo_ref<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options.ref))
  theo_ref <- cbind(iterations, theo_ref)
  final_ref[,2:5] <- final_ref[,2:5]+ theo_ref[,2:5]
  theo_ref['individual'] <- j
  


  options.avg<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=1)
  theo_avg<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options.avg))
  theo_avg <- cbind(iterations, theo_avg)
  final_avg[,2:5] <- final_avg[,2:5]+ theo_avg[,2:5]
  theo_avg['individual'] <- j
  
}


final_ref <- 1/replicate*final_ref
final_avg <- 1/replicate*final_avg


graphConvMC_twokernels(final_ref,final_avg)

graphConvMC_twokernels(final_ref[K1:end,],final_avg[K1:end,])
