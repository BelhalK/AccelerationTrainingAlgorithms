
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/R")
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
  source('mixtureFunctions.R')

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/")
source("/Users/karimimohammedbelhal/Desktop/papers/iem_code/imcem_saemix/plots_se.R")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
###RTTE
timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/rtte_data.csv", header=T, sep=",")
# timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/data/rttellis.csv", header=T, sep=",")
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

K1 = 300
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize50<-50
batchsize25<-25

#Weibull
options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
rtte<-cbind(iterations,rtte)

options_rtteincr<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
rtteincr<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr))
rtteincr<-cbind(iterations,rtteincr)

options_rtteincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
rtteincr25<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr25))
rtteincr25<-cbind(iterations,rtteincr25)

# graphConvMC2_saem(rtte,rtteincr, title="new kernel")


rtte$algo <- 'rwm'
rtteincr$algo <- 'ISAEM50'
rtteincr25$algo <- 'ISAEM25'

rtte_scaled <- rtte[rep(seq_len(nrow(rtte)), each=100/batchsize25),]
rtte_scaled$iterations = 1:(2*(K1+K2+1))

rtte_scaled50 <- rtteincr[rep(seq_len(nrow(rtte)), each=100/batchsize50),]
rtte_scaled50$iterations = 1:(2*(K1+K2+1))

comparison <- 0
# comparison <- rbind(theo_ref,theo_incremental)
comparison <- rbind(rtte_scaled[iterations,],rtte_scaled50[iterations,],rtteincr25)
comparison <- rbind(rtte_scaled[iterations,],rtteincr)

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)

# options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5))
# rttenew<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenew))


replicate = 4
seed0 = 395246

#RWM
final_rwm <- 0
final_50 <- 0
final_25 <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(2,1),c(2,1),c(2,1),c(2,1))
  
  saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(l[[m]],ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

  options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
  rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
  rtte<-cbind(iterations,rtte)
  rtte['individual'] <- m
  rtte$algo <- 'rwm'
  rtte_scaled <- rtte[rep(seq_len(nrow(rtte)), each=100/batchsize25),]
  rtte_scaled$iterations = 1:(2*(K1+K2+1))
  final_rwm <- rbind(final_rwm,rtte)

  options_rtteincr<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
  rtteincr<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr))
  rtteincr<-cbind(iterations,rtteincr)
  rtteincr['individual'] <- m
  rtteincr$algo <- 'ISAEM50'
  rtte_scaled50 <- rtteincr[rep(seq_len(nrow(rtte)), each=100/batchsize50),]
  rtte_scaled50$iterations = 1:(2*(K1+K2+1))
  final_50 <- rbind(final_50,rtteincr)

  options_rtteincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
  rtteincr25<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr25))
  rtteincr25<-cbind(iterations,rtteincr25)
  rtteincr25['individual'] <- m
  rtteincr25$algo <- 'ISAEM25'
  final_25 <- rbind(final_25,rtteincr25)
}

graphConvMC_new(final_rwm, title="RWM")
graphConvMC_diff(final_rwm,final_mix, title="RWM")

#MC STUDY


final_rwm <- 0
final_ref <- 0
var_rwm <- 0
error_rwm <- 0
final_mix <- 0
final_mix25 <- 0
var_mix <- 0
var_mix25 <- 0
error_mix <- 0
error_mix25 <- 0


lambda_true <- 2
o_lambda_true <- 0.3
beta_true <- 2
o_beta_true <- 0.3
final_mix <- 0
true_param <- c(lambda_true,beta_true,o_lambda_true,o_beta_true)
true_param <- data.frame("lambda" = lambda_true, "beta" = beta_true,"omega2.lambda" = o_lambda_true,"omega2.beta" = o_beta_true)

seed0 = 39546
replicate = 20

for (j in 1:replicate){

     model2 <- inlineModel("

  [LONGITUDINAL]
  input = {beta,lambda}  

  EQUATION:
  h=(beta/lambda)*(t/lambda)^(beta-1)

  DEFINITION:
  e = {type               = event, 
       rightCensoringTime = 6,  
       hazard             = h}
  [INDIVIDUAL]
  input={lambda_pop, o_lambda,beta_pop, o_beta}
                        
  DEFINITION:
  lambda  ={distribution=lognormal, prediction=lambda_pop,  sd=o_lambda}
  beta  ={distribution=lognormal, prediction=beta_pop,  sd=o_beta}
       ")


  p <- c(lambda_pop=lambda_true, o_lambda= o_lambda_true,
         beta_pop = beta_true,o_beta = o_beta_true)
  h <- list(name='h', time=seq(0, 6, by=1))
  e <- list(name='e', time=0)

  N <- 50
  res <- simulx(model     = model2, 
                settings  = list(seed=j*123),
                parameter = p, 
                output    = list(h,e), 
                 group     = list(size = N))
  

writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_rtte.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_rtte.csv", header=T, sep=","))
  timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_rtte.csv", header=T, sep=",")
  timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]

  saemix.data<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))


saemix.model<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))
  print(j)

  options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100)
  theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  ML <- theo_ref[,2:5]
  # ML[1:(end+1),]<- theo_ref[end+1,2:5]
  ML[1:(end+1),1:4]<- true_param
  error_rwm <- error_rwm + (theo_ref[,2:5]-ML)^2
  theo_ref['individual'] <- j
  final_ref <- rbind(final_ref,theo_ref)


  options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
  theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
  theo_mix <- cbind(iterations, theo_mix)
  ML <- theo_mix[,2:5]
  # ML[1:(end+1),]<- theo_mix[end+1,2:5]
  ML[1:(end+1),1:4]<- true_param  
  error_mix <- error_mix + (theo_mix[,2:5]-ML)^2
  theo_mix['individual'] <- j
  final_mix <- rbind(final_mix,theo_mix)
  
  options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
  theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
  theo_mix25 <- cbind(iterations, theo_mix25)
  ML <- theo_mix25[,2:5]
  # ML[1:(end+1),]<- theo_mix25[end+1,2:5]
  ML[1:(end+1),1:4]<- true_param 
  error_mix25 <- error_mix25 + (theo_mix25[,2:5]-ML)^2
  theo_mix25['individual'] <- j
  final_mix25 <- rbind(final_mix25,theo_mix25)
}


graphConvMC_diff(final_ref,final_ref,final_ref)
graphConvMC_diff(final_ref,final_mix,final_mix25)

error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix
error_mix25 <- 1/replicate*error_mix25

error_rwm <- cbind(iterations, error_rwm)
error_mix <- cbind(iterations, error_mix)
error_mix25 <- cbind(iterations, error_mix25)

err_mix<- theo_ref
err_rwm<- theo_ref
err_mix25<- theo_ref


err_rwm[,2:5] <- error_rwm[,2:5]
err_mix[,2:5] <- error_mix[,2:5]
err_mix25[,2:5] <- error_mix25[,2:5]


err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = seq(1, 4*end, by=4)

err_mix_scaled <- err_rwm
err_mix_scaled$iterations = seq(1, 2*end, by=2)



# c <- graphConvMC_se2(err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)])
c <- graphConvMC_sec(err_rwm_scaled[2:end,c(1,2,6)],err_mix_scaled[2:end,c(1,2,6)],err_mix25[2:end,c(1,2,6)])
d <- graphConvMC_sed(err_rwm_scaled[2:end,c(1,4,6)],err_mix_scaled[2:end,c(1,4,6)],err_mix25[2:end,c(1,4,6)])

grid.arrange(c,d, ncol=2)


