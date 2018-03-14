
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
source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/plots.R")
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
options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,sampling='seq', map.range=c(0))
rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
rtte<-cbind(iterations,rtte)

options_rtteincr<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,sampling='randompass', map.range=c(0), nb.replacement=batchsize50)
rtteincr50<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr))
rtteincr50<-cbind(iterations,rtteincr50)

options_rtteincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,sampling='randompass', map.range=c(0), nb.replacement=batchsize25)
rtteincr25<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr25))
rtteincr25<-cbind(iterations,rtteincr25)

options_rtteincr75<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,sampling='randompass', map.range=c(0), nb.replacement=75)
rtteincr75<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr75))
rtteincr75<-cbind(iterations,rtteincr75)


options_rtteincr85<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,sampling='randompass', map.range=c(0), nb.replacement=85)
rtteincr85<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr85))
rtteincr85<-cbind(iterations,rtteincr85)
# graphConvMC2_saem(rtte,rtteincr, title="new kernel")

rtte_scaled <- rtte
rtteincr50_scaled <- rtteincr50
rtteincr25_scaled <- rtteincr25
rtteincr75_scaled <- rtteincr75
rtteincr85_scaled <- rtteincr85


rtte_scaled$iterations = rtte_scaled$iterations*1
rtteincr50_scaled$iterations = rtteincr50_scaled$iterations*0.5
rtteincr25_scaled$iterations = rtteincr25_scaled$iterations*0.25
rtteincr75_scaled$iterations = rtteincr75_scaled$iterations*0.75
rtteincr85_scaled$iterations = rtteincr85_scaled$iterations*0.85

graphConvMC_5(rtte_scaled,rtteincr50_scaled,rtteincr25_scaled,rtteincr75_scaled,rtteincr85_scaled)


# options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5))
# rttenew<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenew))


replicate = 3
seed0 = 395246

#RWM
final_rwm <- 0
final_50 <- 0
final_25 <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(2,1),c(1,2),c(3,3),c(2,1))
  
  saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(l[[m]],ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE),omega.init=matrix(c(1/m,0,0,1/m),ncol=2,byrow=TRUE))

  options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
  rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
  rtte<-cbind(iterations,rtte)
  rtte['individual'] <- m
  rtte$algo <- 'rwm'
  rtte_scaled <- rtte[-1,]
  rtte_scaled$iterations = seq(1, 4*end, by=4)
  rtte_scaled <- rtte_scaled[rep(seq_len(nrow(rtte_scaled)), each=4),]
  final_rwm <- rbind(final_rwm,rtte_scaled[0:end,])

  options_rtteincr<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
  rtteincr<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr))
  rtteincr<-cbind(iterations,rtteincr)
  rtteincr['individual'] <- m
  rtteincr$algo <- 'ISAEM50'
  rtte_scaled50 <- rtte[-1,]
  rtte_scaled50$iterations = seq(1, 2*end, by=2)
  rtte_scaled50 <- rtte_scaled50[rep(seq_len(nrow(rtte_scaled50)), each=100/batchsize50),]
  final_50 <- rbind(final_50,rtte_scaled50[0:end,])

  options_rtteincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
  rtteincr25<-data.frame(saemix_incremental(saemix.model_rtte,saemix.data_rtte,options_rtteincr25))
  rtteincr25<-cbind(iterations,rtteincr25)
  rtteincr25['individual'] <- m
  rtteincr25 <- rtteincr25[-1,]
  rtteincr25$algo <- 'ISAEM25'
  rtteincr25$iterations = seq(1, end, by=1)
  rtteincr25[end,] <- rtteincr25[(end-2),]
  final_25 <- rbind(final_25,rtteincr25[0:end,])
}

a <- graphConvMC_diffz(final_rwm[,c(1,2,6)],final_50[,c(1,2,6)],final_25[,c(1,2,6)])
b <- graphConvMC_diffw(final_rwm[,c(1,4,6)],final_50[,c(1,4,6)],final_25[,c(1,4,6)])

grid.arrange(a,b, ncol=2)

#MC STUDY




K1 = 300
K2 = 100
iterations = 0:(K1+K2-1)
end = K1+K2
batchsize50<-50
batchsize25<-25


final_ref <- 0
final_mix50 <- 0
final_mix25 <- 0
final_mix75 <- 0
final_mix85 <- 0

error_rwm <- 0
error_mix50 <- 0
error_mix25 <- 0
error_mix75 <- 0
error_mix85 <- 0


lambda_true <- 2
o_lambda_true <- sqrt(0.3)
beta_true <- 2
o_beta_true <- sqrt(0.3)
final_mix <- 0
true_param <- c(lambda_true,beta_true,o_lambda_true,o_beta_true)
true_param <- data.frame("lambda" = lambda_true, "beta" = beta_true,"omega2.lambda" = o_lambda_true^2,"omega2.beta" = o_beta_true^2)


seed0 = 39546
replicate = 4

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
  h <- list(name='h', time=seq(0, 3, by=1))
  e <- list(name='e', time=0)

  N <- 100
  res <- simulx(model     = model2, 
                parameter = p, 
                output    = list(h,e), 
                 group     = list(size = N))
  
writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_rtte.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_rtte.csv", header=T, sep=","))
  timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_rtte.csv", header=T, sep=",")
  timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
  saemix.data<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))


saemix.model<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(5,5),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))
  print(j)

  options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100, sampling='')
  theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref[-1,])
  ML <- theo_ref[,2:5]
  # ML[1:(end+1),]<- theo_ref[end+1,2:5]
  ML[1:nrow(ML),]<-true_param
  error_rwm <- error_rwm + (theo_ref[,2:5]-ML)^2
  theo_ref['individual'] <- j
  final_ref <- rbind(final_ref,theo_ref)


  options.incremental50<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50, sampling='randompass')
  theo_mix50<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental50))
  theo_mix50 <- cbind(iterations, theo_mix50[-1,])
  
  # ML[1:(end+1),]<- theo_mix50[end+1,2:5]
  error_mix50 <- error_mix50 + (theo_mix50[,2:5]-ML)^2
  theo_mix50['individual'] <- j
  final_mix50 <- rbind(final_mix50,theo_mix50)
  
  options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25, sampling='randompass')
  theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
  theo_mix25 <- cbind(iterations, theo_mix25[-1,])
  
  # ML[1:(end+1),]<- theo_mix25[end+1,2:5]
  error_mix25 <- error_mix25 + (theo_mix25[,2:5]-ML)^2
  theo_mix25['individual'] <- j
  final_mix25 <- rbind(final_mix25,theo_mix25)

  # options.incremental75<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),
  #  nb.replacement=75, sampling='randompass')
  # theo_mix75<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental75))
  # theo_mix75 <- cbind(iterations, theo_mix75[-1,])
  
  # # ML[1:(end+1),]<- theo_mix75[end+1,2:5]
  # error_mix75 <- error_mix75 + (theo_mix75[,2:5]-ML)^2
  # theo_mix75['individual'] <- j
  # final_mix75 <- rbind(final_mix75,theo_mix75)

  # options.incremental85<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),
  #   nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=85, sampling='randompass')
  # theo_mix85<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental85))
  # theo_mix85 <- cbind(iterations, theo_mix85[-1,])
  
  # # ML[1:(end+1),]<- theo_mix85[end+1,2:5]
  # error_mix85 <- error_mix85 + (theo_mix85[,2:5]-ML)^2
  # theo_mix85['individual'] <- j
  # final_mix85 <- rbind(final_mix85,theo_mix85)

}



graphConvMC_diff(final_ref,final_mix75,final_mix85)



final_ref$algo <- '100'
final_mix50$algo <- '50'
final_mix25$algo <- '25'
final_mix75$algo <- '75'
final_mix85$algo <- '85'


final_ref_scaled <- final_ref
final_mix50_scaled <- final_mix50
final_mix25_scaled <- final_mix25
final_mix75_scaled <- final_mix75
final_mix85_scaled <- final_mix85


final_ref_scaled$iterations = final_ref_scaled$iterations*1
final_mix50_scaled$iterations = final_mix50_scaled$iterations*0.5
final_mix25_scaled$iterations = final_mix25_scaled$iterations*0.25
final_mix75_scaled$iterations = final_mix75_scaled$iterations*0.75
final_mix85_scaled$iterations = final_mix85_scaled$iterations*0.85



for (i in 2:5){
comparison <- 0

comparison <- rbind(final_ref_scaled[,c(1,i,6,7)],final_mix25_scaled[,c(1,i,6,7)],
                final_mix50_scaled[,c(1,i,6,7)],final_mix75_scaled[,c(1,i,6,7)],
                final_mix85_scaled[,c(1,i,6,7)])

var <- melt(comparison, id.var = c('iterations','algo','individual'), na.rm = TRUE)

beta0 <- seplot(var,colnames(final_ref_scaled)[i], title="comparison",legend=TRUE)

}


for (i in c(3,7)){
comparison <- 0

comparison <- rbind(final_ref_scaled[,c(1,i,6,7)],
                final_mix85_scaled[,c(1,i,6,7)])

var <- melt(comparison, id.var = c('iterations','algo','individual'), na.rm = TRUE)

beta0 <- seplot(var,colnames(final_ref_scaled)[i], title="comparison",legend=TRUE)

}



error_rwm <- 1/replicate*error_rwm
err_rwm<- theo_ref[,]
err_rwm[,2:5] <- error_rwm[,]
err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = err_rwm_scaled$iterations*1
err_rwm_scaled$algo <- 'SAEM'
err_rwm_scaled$method <- 'seq'


error_mix50 <- 1/replicate*error_mix50
err_mix50<- theo_ref[,]
err_mix50[,2:5] <- error_mix50[,]
err_mix50_scaled <- err_mix50
err_mix50_scaled$iterations = err_mix50_scaled$iterations*0.5
err_mix50_scaled$algo <- 'ISAEM 50'
err_mix50_scaled$method <- 'randompass'

error_mix25 <- 1/replicate*error_mix25
err_mix25<- theo_ref[,]
err_mix25[,2:5] <- error_mix25[,]
err_mix25_scaled <- err_mix25
err_mix25_scaled$iterations = err_mix25_scaled$iterations*0.25
err_mix25_scaled$algo <- 'ISAEM 25'
err_mix25_scaled$method <- 'randompass'

error_mix75 <- 1/replicate*error_mix75
err_mix75<- theo_ref[,]
err_mix75[,2:5] <- error_mix75[,]
err_mix75_scaled <- err_mix75
err_mix75_scaled$iterations = err_mix75_scaled$iterations*0.75
err_mix75_scaled$algo <- 'ISAEM 75'
err_mix75_scaled$method <- 'randompass'

error_mix85 <- 1/replicate*error_mix85
err_mix85<- theo_ref[,]
err_mix85[,2:5] <- error_mix85[,]
err_mix85_scaled <- err_mix85
err_mix85_scaled$iterations = err_mix85_scaled$iterations*0.85
err_mix85_scaled$algo <- 'ISAEM 85'
err_mix85_scaled$method <- 'randompass'





for (i in 2:5){
# i = 6
comparison <- 0
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,7,8)],
                    err_mix50_scaled[0:end,c(1,i,7,8)],
                    err_mix25_scaled[0:end,c(1,i,7,8)])

var <- melt(comparison, id.var = c('iterations','algo','method'), na.rm = TRUE)


prec <- seplot2(var,colnames(err_rwm_scaled)[i], title="ALGO - EM (same complexity)",legend=TRUE)
# setwd("/Users/karimimohammedbelhal/Desktop/")
# ggsave(paste("precwarfa_", i, ".png", sep=""),prec)
}

