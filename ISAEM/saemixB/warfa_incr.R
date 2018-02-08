
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

  source('/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/R/mixtureFunctions.R')
  source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/plots.R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB")


library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)


model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))



K1 = 400
K2 = 100
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

ka_true <- 0.7272239
V_true <- 7.536013
k_true <- 1.756585
o_ka <- sqrt(0.5078622)
o_V <- sqrt(0.03980655)
o_k <- sqrt(0.04624554)
a_true<-1

true_param <- data.frame("ka" = ka_true, "V" = V_true, "k" = k_true, "omega2.ka"=o_ka ,"omega2.V"= o_V,"omega2.k"= o_k, "a" = a_true)
ML<- true_param[rep(seq_len(nrow(true_param)), (end)), ]
seed0 = 39546
replicate = 20


for (m in 1:replicate){
  print(m)
model2 <- inlineModel("

                   [INDIVIDUAL]
                  input = {ka_pop, omega_ka, V_pop, omega_V, k_pop, omega_k}

                  DEFINITION:
                  ka = {distribution=lognormal, typical=ka_pop, sd=omega_ka}
                  V = {distribution=lognormal, typical=V_pop, sd=omega_V}
                  k = {distribution=lognormal, typical=k_pop, sd=omega_k} 
                  [LONGITUDINAL]
                  input = {ka,V, k, a}
                  EQUATION:
                  Cc = pkmodel(ka,V,k)
                  DEFINITION:
                  y = {distribution=normal, prediction=Cc, sd=a}
                          ")

    adm <- list(time=0, amount=100)
    y <- list(name='y', time=seq(1, 10, by=1))
    p <- c(ka_pop=ka_true, omega_ka=o_ka,
           V_pop=V_true, omega_V=o_V, 
           k_pop=k_true, omega_k=o_k,
           a=a_true)
    g <- list(size=100, level="individual")
    res<-simulx(model=model2, parameter=p, output=y, treatment=adm, group=g)
  
  warfarin.saemix <- res$y
  warfarin.saemix$amount <- 100

  saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(5,15,5,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(2,0,0,0,2,0,0,0,2),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


  saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")


  options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='')
  theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
  # ML[1:(end+1),]<- theo_ref[end+1,2:8]
  error_rwm <- error_rwm + (theo_ref[-1,]-ML)^2
  final_ref <- rbind(final_ref,theo_ref)

  #SEQ
  options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='seq')
  theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
  # ML[1:(end+1),]<- theo_mix[end+1,2:8]
  error_mixseq <- error_mixseq + (theo_mix[-1,]-ML)^2
  final_mix <- rbind(final_mix,theo_mix)

  options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='seq')
  theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
  # ML[1:(end+1),]<- theo_mix25[end+1,2:8]
  error_mix25seq <- error_mix25seq + (theo_mix25[-1,]-ML)^2
  final_mix25 <- rbind(final_mix25,theo_mix25)

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

  #RANDOMITER
  options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randomiter')
  theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
  # ML[1:(end+1),]<- theo_mix[end+1,2:8]
  error_mixiter <- error_mixiter + (theo_mix[-1,]-ML)^2
  final_mix <- rbind(final_mix,theo_mix)

  options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randomiter')
  theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
  # ML[1:(end+1),]<- theo_mix25[end+1,2:8]
  error_mix25iter <- error_mix25iter + (theo_mix25[-1,]-ML)^2
  final_mix25 <- rbind(final_mix25,theo_mix25)
 
}

# graphConvMC_diff(final_ref,final_ref,final_ref)
# graphConvMC_diff(final_mix25,final_mix25,final_mix25)
# graphConvMC_diff(final_ref,final_mix,final_mix25)

error_rwm <- 1/replicate*error_rwm
error_mixseq <- 1/replicate*error_mixseq
error_mix25seq <- 1/replicate*error_mix25seq

error_mixpass <- 1/replicate*error_mixpass
error_mix25pass <- 1/replicate*error_mix25pass

error_mixiter <- 1/replicate*error_mixiter
error_mix25iter <- 1/replicate*error_mix25iter


error_rwm <- cbind(iterations,error_rwm)
error_mixseq <- cbind(iterations,error_mixseq)
error_mix25seq <- cbind(iterations,error_mix25seq)

error_mixpass <- cbind(iterations,error_mixpass)
error_mix25pass <- cbind(iterations,error_mix25pass)

error_mixiter <- cbind(iterations,error_mixiter)
error_mix25iter <- cbind(iterations,error_mix25iter)


err_rwm_scaled <- data.frame(error_rwm)
err_rwm_scaled$iterations = seq(1, 4*end, by=4)

err_mixseq_scaled <- data.frame(error_mixseq)
err_mixseq_scaled$iterations = seq(1, 2*end, by=2)

err_mixpass_scaled <- data.frame(error_mixpass)
err_mixpass_scaled$iterations = seq(1, 2*end, by=2)

err_mixiter_scaled <- data.frame(error_mixiter)
err_mixiter_scaled$iterations = seq(1, 2*end, by=2)


error_mix25seq <- data.frame(error_mix25seq)
error_mix25pass <- data.frame(error_mix25pass)
error_mix25iter <- data.frame(error_mix25iter)

error_mix25seq$iterations = 1:((K1+K2))
error_mix25pass$iterations = 1:((K1+K2))
error_mix25iter$iterations = 1:((K1+K2))




err_rwm_scaled$algo <- 'SAEM'
err_rwm_scaled$method <- 'seq'

err_mixseq_scaled$algo <- 'ISAEM50'
err_mixpass_scaled$algo <- 'ISAEM50'
err_mixiter_scaled$algo <- 'ISAEM50'

err_mixseq_scaled$method <- 'seq'
err_mixpass_scaled$method <- 'pass'
err_mixiter_scaled$method <- 'iter'


error_mix25seq$algo <- 'ISAEM25'
error_mix25pass$algo <- 'ISAEM25'
error_mix25iter$algo <- 'ISAEM25'

error_mix25seq$method <- 'seq'
error_mix25pass$method <- 'pass'
error_mix25iter$method <- 'iter'



for (i in 2:8){
# i = 6
comparison <- 0
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,9,10)],err_mixseq_scaled [0:end,c(1,i,9,10)],error_mix25seq[0:end,c(1,i,9,10)],
                                              err_mixpass_scaled [0:end,c(1,i,9,10)],error_mix25pass[0:end,c(1,i,9,10)],
                                              err_mixiter_scaled [0:end,c(1,i,9,10)],error_mix25iter[0:end,c(1,i,9,10)])

var <- melt(comparison, id.var = c('iterations','algo','method'), na.rm = TRUE)


prec <- seplot(var, title="ALGO - EM (same complexity)",legend=TRUE)
# setwd("/Users/karimimohammedbelhal/Desktop/")
# ggsave(paste("precwarfa_", i, ".png", sep=""),prec)
}


# c <- graphConvMC_se2(err_rwm_scaled[-1,c(1,2,8)],err_mix_scaled[-1,c(1,2,8)],err_mix25[-1,c(1,2,8)])
# d <- graphConvMC_se2(err_rwm_scaled[-1,c(1,3,8)],err_mix_scaled[-1,c(1,3,8)],err_mix25[-1,c(1,3,8)])

# grid.arrange(c,d, ncol=2)

# c <- graphConvMC_se2(err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)])
c <- graphConvMC_sec(err_rwm_scaled[0:end,c(1,2,9)],err_mix_scaled [0:end,c(1,2,9)],err_mix25[0:end,c(1,2,9)])
d <- graphConvMC_sed(err_rwm_scaled[0:end,c(1,5,9)],err_mix_scaled[0:end,c(1,5,9)],err_mix25[0:end,c(1,5,9)])

grid.arrange(c,d, ncol=2)

e <- graphConvMC_sec(err_rwm_scaled[0:end,c(1,3,9)],err_mix_scaled[0:end,c(1,3,9)],err_mix25[0:end,c(1,3,9)])
f <- graphConvMC_sed(err_rwm_scaled[0:end,c(1,6,9)],err_mix_scaled[0:end,c(1,6,9)],err_mix25[0:end,c(1,6,9)])

grid.arrange(e,f, ncol=2)

g <- graphConvMC_sec(err_rwm_scaled[0:end,c(1,4,9)],err_mix_scaled[0:end,c(1,4,9)],err_mix25[0:end,c(1,4,9)])
h <- graphConvMC_sed(err_rwm_scaled[0:end,c(1,7,9)],err_mix_scaled[0:end,c(1,7,9)],err_mix25[0:end,c(1,7,9)])

grid.arrange(g,h, ncol=2)


