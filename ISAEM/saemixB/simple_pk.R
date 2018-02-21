
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



warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  k <- Cl/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*time)-exp(-ka*time))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


K1 = 1000
K2 = 300
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
  l = list(c(1,7,1,0,0,0),c(0.8,7.2,0.8,0,0,0),c(1.2,6.8,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

  options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='seq')
  warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
  warfa<-cbind(iterations,warfa)
  warfa['individual'] <- m
  warfa$algo <- 'rwm'
  warfa_scaled <- warfa[-1,]
  warfa_scaled$iterations = seq(1, 4*end, by=4)
  warfa_scaled <- warfa_scaled[rep(seq_len(nrow(warfa_scaled)), each=100/batchsize25),]
  final_rwm <- rbind(final_rwm,warfa_scaled[0:end,])

  options_warfaincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randompass')
  warfaincr25<-data.frame(saemix_incremental(saemix.model_warfa,saemix.data_warfa,options_warfaincr25))
  warfaincr25<-cbind(iterations,warfaincr25)
  warfaincr25['individual'] <- m
  warfaincr25$algo <- 'ISAEM25'
  warfaincr25 <- warfaincr25[-1,]
  warfaincr25$iterations = seq(1, end, by=1)
  final_25 <- rbind(final_25,warfaincr25[0:end,])

  options_warfaincr50<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randompass')
  warfaincr50<-data.frame(saemix_incremental(saemix.model_warfa,saemix.data_warfa,options_warfaincr50))
  warfaincr50<-cbind(iterations,warfaincr50)
  warfaincr50['individual'] <- m
  warfaincr50$algo <- 'ISAEM50'
  warfa_scaled50 <- warfaincr50[-1,]
  warfa_scaled50$iterations = seq(1, 2*end, by=2)
  warfa_scaled50 <- warfa_scaled50[rep(seq_len(nrow(warfa_scaled50)), each=100/batchsize50),]
  final_50 <- rbind(final_50,warfa_scaled50[0:end,])
}

graphConvMC2_saem(warfa[-1,c(1:8)],warfa[-1,c(1:8)], title="new kernel") 
graphConvMC2_saem(warfa[-1,c(1:8)],warfaincr25[,c(1:8)], title="new kernel") 

warfa$algo <- 'rwm'
warfaincr25$algo <- 'ISAEM25'
warfaincr50$algo <- 'ISAEM50'

warfa_scaled <- warfa[rep(seq_len(nrow(warfa)), each=100/batchsize25),]
warfa_scaled$iterations = 1:(4*(K1+K2+1))

warfa_scaled50 <- warfaincr50[rep(seq_len(nrow(warfaincr50)), each=100/batchsize50),]
warfa_scaled50$iterations = 1:(2*(K1+K2+1))

comparison <- 0
# comparison <- rbind(theo_ref,theo_incremental)
# comparison <- rbind(warfa_scaled[iterations,],warfaincr)
comparison <- rbind(warfa_scaled[iterations,],warfa_scaled50[iterations,],warfaincr25)

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)

K1 = 600
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

ka_true <- 2
V_true <- 13
Cl_true <- 1
o_ka <- sqrt(0.2)
o_V <- sqrt(0.1)
o_Cl <- sqrt(0.3)
a_true<-1


true_param <- data.frame("ka" = ka_true, "V" = V_true, "Cl" = Cl_true, "omega2.ka"=o_ka^2 ,"omega2.V"= o_V^2,"omega2.Cl"= o_Cl^2, "a" = a_true)
ML<- true_param[rep(seq_len(nrow(true_param)), (end)), ]
seed0 = 39546
replicate = 3


for (m in 1:replicate){
  print(m)

myModel <- inlineModel("


[INDIVIDUAL]
input = {ka_pop, V_pop, Cl_pop, omega_ka, omega_V, omega_Cl}
DEFINITION:
ka = {distribution=lognormal, reference=ka_pop, sd=omega_ka}
V  = {distribution=lognormal, reference=V_pop,  sd=omega_V }
Cl = {distribution=lognormal, reference=Cl_pop, sd=omega_Cl}


[LONGITUDINAL]
input = {ka, V, Cl,a}
EQUATION:
C = pkmodel(ka,V,Cl)
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=400
pop.param   <- c(
  ka_pop  = ka_true,    omega_ka  = o_ka,
  V_pop   = V_true,   omega_V   = o_V,
  Cl_pop  = Cl_true,    omega_Cl  = o_Cl, a =a_true)
  
res <- simulx(model     = myModel,
              parameter = pop.param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(0,2,by=1)))
  
  # writeDatamlx(res, result.file = "res.csv")
  # head(read.csv("res.csv"))
  # tab <- read.csv("res.csv")
  
  warfarin.saemix <- res$y
  warfarin.saemix$amount <- 100

  
  saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(5,10,5,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

  saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")


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
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,9,10)],err_mixpass_scaled [0:end,c(1,i,9,10)],error_mix25pass[0:end,c(1,i,9,10)])

var <- melt(comparison, id.var = c('iterations','algo','method'), na.rm = TRUE)


prec <- seplot(var, title="ALGO - EM (same complexity)",legend=TRUE)
# setwd("/Users/karimimohammedbelhal/Desktop/")
# ggsave(paste("precwarfa_", i, ".png", sep=""),prec)
}

