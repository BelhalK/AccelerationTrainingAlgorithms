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

library("rJava")
library("rCMA")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

# Doc
# Doc
data(cow.saemix)
saemix.data<-saemixData(name.data=cow.saemix,header=TRUE,name.group=c("cow"), 
  name.predictors=c("time"),name.response=c("weight"), 
  name.covariates=c("birthyear","twin","birthrank"), 
  units=list(x="days",y="kg",covariates=c("yr","-","-")))


# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
growthcow<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, a, b, k)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  a<-psi[id,1]
  b<-psi[id,2]
  k<-psi[id,3]
  f<-a*(1-b*exp(-k*x))
  return(f)
}
saemix.model<-saemixModel(model=growthcow,type="structural",
  description="Exponential growth model", 
  psi0=matrix(c(700,0.9,0.02,0,0,0),ncol=3,byrow=TRUE, 
  dimnames=list(NULL,c("A","B","k"))),transform.par=c(1,1,1),fixed.estim=c(1,1,1), 
  covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
  omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")



K1 = 300
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50
replicate = 3
seed0 = 395246

#RWM
final_rwm <- 0
final_50 <- 0
final_25 <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(500,3,0.02),c(600,2,0.02),c(700,4,0.02),c(700,0.9,0.02))


saemix.model<-saemixModel(model=growthcow,type="structural",
  description="Exponential growth model", 
  psi0=matrix(l[[m]],ncol=3,byrow=TRUE, 
  dimnames=list(NULL,c("A","B","k"))),transform.par=c(1,1,1),fixed.estim=c(1,1,1), 
  covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
  omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE),error.model="constant")

  options_cow<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
  cow<-data.frame(saemix(saemix.model,saemix.data,options_cow))
  cow<-cbind(iterations,cow)
  cow['individual'] <- m
  cow$algo <- 'rwm'
  cow_scaled <- cow[-1,]
  cow_scaled$iterations = seq(1, 4*end, by=4)
  cow_scaled <- cow_scaled[rep(seq_len(nrow(cow_scaled)), each=4),]
  final_rwm <- rbind(final_rwm,cow_scaled[0:end,])

  options_cowincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
  cowincr25<-data.frame(saemix_incremental(saemix.model,saemix.data,options_cowincr25))
  cowincr25<-cbind(iterations,cowincr25)
  cowincr25['individual'] <- m
  cowincr25$algo <- 'ISAEM25'
  cowincr25 <- cowincr25[-1,]
  cowincr25$iterations = seq(1, end, by=1)
  final_25 <- rbind(final_25,cowincr25[0:end,])

  options_cowincr50<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
  cowincr50<-data.frame(saemix_incremental(saemix.model,saemix.data,options_cowincr50))
  cowincr50<-cbind(iterations,cowincr50)
  cowincr50['individual'] <- m
  cowincr50$algo <- 'ISAEM50'
  cow_scaled50 <- cowincr50[-1,]
  cow_scaled50$iterations = seq(1, 2*end, by=2)
  cow_scaled50 <- cow_scaled50[rep(seq_len(nrow(cow_scaled50)), each=2),]
  final_50 <- rbind(final_50,cow_scaled50[0:end,])
}


a <- graphConvMC_diffz(final_rwm[,c(1,2,9)],final_50[,c(1,2,9)],final_25[,c(1,2,9)])
b <- graphConvMC_diffw(final_rwm[,c(1,6,9)],final_50[,c(1,6,9)],final_25[,c(1,6,9)])

grid.arrange(a,b, ncol=2)




K1 = 500
K2 = 300
iterations = 1:(K1+K2)
end = K1+K2
batchsize25 = 25
batchsize50 = 50


final_rwm <- 0
final_ref <- 0
final_mix <- 0
final_mix25 <- 0

error_rwm <- 0
error_mixseq <- 0
error_mix25seq <- 0

error_mixiter <- 0
error_mix25iter <- 0

error_mixpass <- 0
error_mix25pass <- 0

# ka_true <- 0.7272239
# V_true <- 7.536013
# k_true <- 1.756585
# o_ka <- sqrt(0.5078622)
# o_V <- sqrt(0.03980655)
# o_k <- sqrt(0.04624554)
# a_true<-1

d_true <- 200
b_true <- 3
k_true <- 0.5
o_d <- sqrt(0.3)
o_b <- sqrt(0.2)
o_k <- sqrt(0.2)
a_true<-0.1


true_param <- data.frame("d" = d_true, "b" = b_true, "k" = k_true, "omega2.d"=o_d^2 ,"omega2.b"= o_b^2,"omega2.k"= o_k^2, "a" = a_true)
ML<- true_param[rep(seq_len(nrow(true_param)), (end)), ]
seed0 = 39546
replicate = 6

growthcow<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, a, b, k)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  d<-psi[id,1]
  b<-psi[id,2]
  k<-psi[id,3]
  f<-d*(1-b*exp(-k*x))
  return(f)
}

for (m in 1:replicate){
  print(m)

myModel <- inlineModel("


[INDIVIDUAL]
input = {d_pop, b_pop, k_pop, omega_d, omega_b, omega_k}
DEFINITION:
d = {distribution=lognormal, reference=d_pop, sd=omega_d}
b  = {distribution=lognormal, reference=b_pop,  sd=omega_b}
k = {distribution=lognormal, reference=k_pop, sd=omega_k}


[LONGITUDINAL]
input = {d, b, k,a}
EQUATION:
C = d*(1-b*exp(-k*t))
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=200
pop.param   <- c(
  d_pop  = d_true,    omega_d  = o_d,
  b_pop   = b_true,   omega_b   = o_b,
  k_pop  = k_true,    omega_k  = o_k, a =a_true)
  
res <- simulx(model     = myModel,
              parameter = pop.param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(0,5,by=1)))
  
  # writeDatamlx(res, result.file = "res.csv")
  # head(read.csv("res.csv"))
  # tab <- read.csv("res.csv")
  
  cow.saemix <- res$y

  head(cow.saemix)

  saemix.model<-saemixModel(model=growthcow,type="structural",
  description="Exponential growth model", 
  psi0=matrix(c(130,0.9,0.02,0,0,0),ncol=3,byrow=TRUE, 
  dimnames=list(NULL,c("d","b","k"))),transform.par=c(1,1,1),fixed.estim=c(1,1,1), 
  covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
  omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")

saemix.data<-saemixData(name.data=cow.saemix,header=TRUE,name.group=c("id"), 
  name.predictors=c("time"),name.response=c("y"))



  options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='seq')
  theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
  # ML[1:(end+1),]<- theo_ref[end+1,2:8]
  error_rwm <- error_rwm + (theo_ref[-1,]-ML)^2
  final_ref <- rbind(final_ref,theo_ref)

  # #SEQ
  # options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='seq')
  # theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
  # # ML[1:(end+1),]<- theo_mix[end+1,2:8]
  # error_mixseq <- error_mixseq + (theo_mix[-1,]-ML)^2
  # final_mix <- rbind(final_mix,theo_mix)

  # options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='seq')
  # theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
  # # ML[1:(end+1),]<- theo_mix25[end+1,2:8]
  # error_mix25seq <- error_mix25seq + (theo_mix25[-1,]-ML)^2
  # final_mix25 <- rbind(final_mix25,theo_mix25)

  #RANDOMPASS
  options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randompass')
  theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
  # ML[1:(end+1),]<- theo_mix[end+1,2:8]
  error_mixpass <- error_mixpass + (theo_mix[-1,]-ML)^2
  final_mix <- rbind(final_mix,theo_mix)

  options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randompass')
  theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
  # ML[1:(end+1),]<- theo_mix25[end+1,2:8]
  error_mix25pass <- error_mix25pass + (theo_mix25[-1,]-ML)^2
  final_mix25 <- rbind(final_mix25,theo_mix25)

  # #RANDOMITER
  # options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randomiter')
  # theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
  # # ML[1:(end+1),]<- theo_mix[end+1,2:8]
  # error_mixiter <- error_mixiter + (theo_mix[-1,]-ML)^2
  # final_mix <- rbind(final_mix,theo_mix)

  # options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randomiter')
  # theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
  # # ML[1:(end+1),]<- theo_mix25[end+1,2:8]
  # error_mix25iter <- error_mix25iter + (theo_mix25[-1,]-ML)^2
  # final_mix25 <- rbind(final_mix25,theo_mix25)
 
}




error_rwm <- 1/replicate*error_rwm

error_mixpass <- 1/replicate*error_mixpass
error_mix25pass <- 1/replicate*error_mix25pass

error_rwm <- cbind(iterations,error_rwm)
error_mixpass <- cbind(iterations,error_mixpass)
error_mix25pass <- cbind(iterations,error_mix25pass)

err_rwm_scaled <- data.frame(error_rwm)
err_rwm_scaled$iterations = seq(1, 4*end, by=4)

err_rwm_scaled$algo <- 'SAEM'
err_rwm_scaled$method <- 'seq'

error_mix25pass <- data.frame(error_mix25pass)
error_mix25pass$iterations = 1:((K1+K2))
error_mix25pass$algo <- 'ISAEM25'
error_mix25pass$method <- 'pass'


err_mixpass_scaled <- data.frame(error_mixpass)
err_mixpass_scaled$iterations = seq(1, 2*end, by=2)
err_mixpass_scaled$algo <- 'ISAEM50'
err_mixpass_scaled$method <- 'pass'

for (i in 2:8){
# i = 6
comparison <- 0
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,9,10)],err_mixpass_scaled [0:end,c(1,i,9,10)],error_mix25pass[0:end,c(1,i,9,10)],
                                              err_mixpass_scaled [0:end,c(1,i,9,10)],error_mix25pass[0:end,c(1,i,9,10)],
                                              err_mixpass_scaled [0:end,c(1,i,9,10)],error_mix25pass[0:end,c(1,i,9,10)])

var <- melt(comparison, id.var = c('iterations','algo','method'), na.rm = TRUE)


prec <- seplot(var, title="ALGO - EM (same complexity)",legend=TRUE)
# setwd("/Users/karimimohammedbelhal/Desktop/")
# ggsave(paste("precwarfa_", i, ".png", sep=""),prec)
}

