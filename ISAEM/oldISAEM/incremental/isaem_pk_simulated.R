# load("warfa_isaem_newkernel.RData")
# save.image("warfa_isaem_newkernel.RData")
# save.image("warfa_isaem_LARGEN.RData")

####20 CHAINS
load("warfa_isaem.RData")
# save.image("warfa_isaem.RData")

# load("warfa_onlinesaem.RData")
# save.image("warfa_onlinesaem.RData")

# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/incremental/R")
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
  
  source('main_incremental.R')
  source('main_estep_incremental.R')


  source('/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/incremental/R/mixtureFunctions.R')
  source("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/incremental/plots.R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/incremental")


library("mlxR")
library("rlist")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)


# warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/csda_new/data/warfarin_data.txt", header=T)



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

model <- inlineModel("


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

N=100

param   <- c(
  ka_pop  = 2,    omega_ka  = 0.3,
  V_pop   = 10,   omega_V   = 0.2,
  Cl_pop  = 1,    omega_Cl  = 0.3, a =1)
  
res <- simulx(model     = model,
              parameter = param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(1,12,by=1)))

 warfarin.saemix <- res$y
 warfarin.saemix$amount <- 100
 saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")

# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(3,3,0.1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),fixed.estim=c(1,1,1),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


# warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/incremental/data/warfarin_data.txt", header=T)
# saemix.data<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")


# saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(3,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))



K1 = 200
K2 = 30
iterations = 0:(K1+K2-1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50 

seed0=3456

nchains = 1
gamma = 1
options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains,nbiter.mcmc = c(2,2,2,0), 
  nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, 
  map.range=c(0), nb.replacement=100,sampling='randomiter',gamma=gamma, algo="full")
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])


# options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, 
#   nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),
#   nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randomiter',gamma=gamma)
# theo_mix75<-saemix_incremental(saemix.model,saemix.data,options.incremental75)
# theo_mix75 <- data.frame(theo_mix75$param)
# theo_mix75 <- cbind(iterations, theo_mix75[-1,])


# options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,0), 
#                           nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,
#                           nbiter.burn =0, nb.replacement=50,sampling='seq',gamma=gamma)
# theo_mix50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
# theo_mix50 <- data.frame(theo_mix50$param)
# theo_mix50 <- cbind(iterations, theo_mix50[-1,])

##############ONLINE SAEM################
options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=25,sampling='randompass',gamma=gamma,algo="online")
theo_mix25online<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25online <- data.frame(theo_mix25online$param)
theo_mix25online <- cbind(iterations, theo_mix25online[-1,])

theo_mix25online_scaled <- theo_mix25online
theo_mix25online_scaled$iterations = theo_mix25online_scaled$iterations*0.25


##############VR SAEM################
options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=25,sampling='randompass',gamma=gamma,algo="vr")
theo_mix25vr<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25vr <- data.frame(theo_mix25vr$param)
theo_mix25vr <- cbind(iterations, theo_mix25vr[-1,])

theo_mix25vr_scaled <- theo_mix25vr
theo_mix25vr_scaled$iterations = theo_mix25vr_scaled$iterations*0.25




options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='seq',gamma=gamma,algo="minibatch")
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])



theo_ref_scaled <- theo_ref
theo_mix25_scaled <- theo_mix25
theo_mix50_scaled <- theo_mix50
theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5


theo_mix75_scaled <- theo_mix75
theo_mix75_scaled$iterations = theo_mix75_scaled$iterations*0.75
theo_mix75_scaled_20chains <- theo_mix75_scaled




# theo_ref_scaled_20chains <- theo_ref_scaled
# theo_mix25_scaled_20chains <- theo_mix25_scaled
# theo_mix50_scaled_20chains <- theo_mix50_scaled

# theo_ref_scaled_10chains <- theo_ref_scaled
# theo_mix25_scaled_10chains <- theo_mix25_scaled
# theo_mix50_scaled_10chains <- theo_mix50_scaled
graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix25_scaled)
graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix25online_scaled)
graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix25vr_scaled)

graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix50online_scaled)

graphConvMC_threekernels(theo_ref_scaled_20chains,theo_mix25_scaled_20chains,theo_mix50_scaled_20chains)
graphConvMC_threekernels(theo_ref_scaled_10chains,theo_mix25_scaled_10chains,theo_mix50_scaled_10chains)
graphConvMC_threekernels(theo_ref_scaled_1chain,theo_mix25_scaled_1chain,theo_mix50_scaled_1chain)


graphConvMC_5(theo_ref_scaled_20chains,theo_mix25_scaled_20chains,theo_mix50_scaled_20chains,theo_mix75_scaled_20chains,theo_mix75_scaled_20chains)
graphConvMC_5(theo_ref_scaled_10chains,theo_mix25_scaled_10chains,theo_mix50_scaled_10chains,theo_mix75_scaled_10chains,theo_mix75_scaled_10chains)
graphConvMC_5(theo_ref_scaled_1chain,theo_mix25_scaled_1chain,theo_mix50_scaled_1chain,theo_mix75_scaled_1chain,theo_mix75_scaled_1chain)
# graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled)
# graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix50_scaled,theo_mix75_scaled)

###NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL####
nchains = 1
gamma = 1
options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains,
  nbiter.mcmc = c(2,2,2,2),nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,
  nbiter.burn = 0, 
 map.range=c(1:4), nb.replacement=100,sampling='seq', gamma = gamma)
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])



options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,2), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:10),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass',gamma=gamma)
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, 
  nbiter.mcmc = c(2,2,2,2), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:4),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass',gamma=gamma)
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])


options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,2), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:10),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=75,sampling='randompass',gamma=gamma)
theo75<-saemix_incremental(saemix.model,saemix.data,options.incremental75)
theo_mix75 <- data.frame(theo75$param)
theo_mix75 <- cbind(iterations, theo_mix75[-1,])



theo_ref_scaled <- theo_ref
theo_mix50_scaled <- theo_mix50
theo_mix75_scaled <- theo_mix75
theo_mix25_scaled <- theo_mix25
theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix75_scaled$iterations = theo_mix75_scaled$iterations*0.75
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25

# graphConvMC_threekernels(theo_ref_scaled,theo_mix50_scaled,theo_mix75_scaled)
graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix75_scaled,theo_mix75_scaled)
# graphConvMC_threekernels(theo_ref,theo_mix50,theo_mix75)
###NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL####


options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
theo_mix75<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental75))
theo_mix75 <- cbind(iterations, theo_mix75[-1,])


options.incremental85<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=85,sampling='randompass')
theo_mix85<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental85))
theo_mix85 <- cbind(iterations, theo_mix85[-1,])


theo_ref_scaled <- theo_ref
theo_mix50_scaled <- theo_mix50
theo_mix25_scaled <- theo_mix25
theo_mix75_scaled <- theo_mix75
theo_mix85_scaled <- theo_mix85


theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25
theo_mix75_scaled$iterations = theo_mix75_scaled$iterations*0.75
theo_mix85_scaled$iterations = theo_mix85_scaled$iterations*0.85

# graphConvMC_threekernels(theo_ref_scaled,theo_mix50_scaled,theo_mix50_scaled)
# graphConvMC_threekernels(theo_ref_scaled,theo_mix50_scaled,theo_mix25_scaled)
graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix75_scaled,theo_mix85_scaled)



final_ref <- 0
final_mix <- 0

ka_true <- 1
V_true <- 10
Cl_true <- 1
a_true <- 1

o_ka <- 0.3 #o^2=0.09
o_V <- 0.2  #o^2=0.04
o_Cl <- 0.3  #o^2=0.09

error_rwm <- 0
error_mix <- 0
error_mix2 <- 0
error_mix20 <- 0


K1 = 500
K2 = 100
iterations = 0:(K1+K2-1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50



true_param <- data.frame("ka" = ka_true, "V" = V_true, "Cl" = Cl_true, "omega2.ka"=o_ka^2 ,"omega2.V"= o_V^2,"omega2.Cl"= o_Cl^2, "a" = a_true)
seed0 = 39546
replicate = 20
for (m in 1:replicate){

model <- inlineModel("


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

N=100

param   <- c(
  ka_pop  = ka_true,    omega_ka  = o_ka,
  V_pop   = V_true,   omega_V   = o_V,
  Cl_pop  = Cl_true,    omega_Cl  = o_Cl, a =1)
  
res <- simulx(model     = model,
              parameter = param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(1,10,by=1)))

 warfarin.saemix <- res$y
 warfarin.saemix$amount <- 100
 saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")

# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,3,0.1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),fixed.estim=c(0,1,0),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,2),
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, 
 map.range=c(1:4), nb.replacement=100,sampling='seq', gamma= 0)
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])
ML <- theo_ref[,2:8]
ML[0:end,1:7]<- theo_ref[end,2:8]
# ML[0:end,1:7]<- true_param
error_rwm <- error_rwm + (theo_ref[,2:8]-ML)^2
theo_ref['individual'] <- m
final_ref <- rbind(final_ref,theo_ref)


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,2), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:4),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass', gamma = 0.2)
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])
ML <- theo_ref[,2:8]
ML[0:end,1:7]<- theo_mix50[end,2:8]
# ML[0:end,1:7]<- true_param
error_mix <- error_mix + (theo_mix50[,2:8]-ML)^2
theo_mix50['individual'] <- m
final_mix <- rbind(final_mix,theo_mix50)


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,2), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:4),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass', gamma = 2)
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])
ML <- theo_ref[,2:8]
ML[0:end,1:7]<- theo_mix50[end,2:8]
# ML[0:end,1:7]<- true_param
error_mix2 <- error_mix2 + (theo_mix50[,2:8]-ML)^2
theo_mix50['individual'] <- m
final_mix <- rbind(final_mix,theo_mix50)

options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,2), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:4),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass', gamma = 20)
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])
ML <- theo_ref[,2:8]
ML[0:end,1:7]<- theo_mix50[end,2:8]
# ML[0:end,1:7]<- true_param
error_mix20 <- error_mix20 + (theo_mix50[,2:8]-ML)^2
theo_mix50['individual'] <- m
final_mix <- rbind(final_mix,theo_mix50)

}



error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix

err_rwm<- theo_ref[-1,]
err_mix<- theo_ref[-1,]

err_rwm[,2:8] <- error_rwm[-1,]
err_mix[,2:8] <- error_mix[-1,]



err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = seq(1, 4*(end-1), by=4)

err_mix_scaled <- err_mix
err_mix_scaled$iterations = seq(1, 2*(end-1), by=2)

err_rwm_scaled$algo <- 'SAEM'
err_rwm_scaled$method <- 'randiter'

err_mix_scaled$algo <- 'ISAEM50 gamma 0.2'
err_mix_scaled$method <- 'randiter'



error_mix50_gamma2 <- 1/replicate*error_mix2
err_mix50_gamma2<- theo_ref[-1,]
err_mix50_gamma2[,2:8] <- error_mix50_gamma2[-1,]
err_mix50_gamma2$iterations = seq(1, (end-1), by=1)
err_mix50_gamma2$algo <- 'ISAEM50 gamma 2'
err_mix50_gamma2$method <- 'randiter'


error_mix50_gamma20 <- 1/replicate*error_mix20
err_mix50_gamma20<- theo_ref[-1,]
err_mix50_gamma20[,2:8] <- error_mix50_gamma20[-1,]
err_mix50_gamma20$iterations = seq(1, (end-1), by=1)
err_mix50_gamma20$algo <- 'ISAEM50 gamma 20'
err_mix50_gamma20$method <- 'randiter'

graphConvMC_se1 <- function(df,df2, df3,df4,title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  df3$individual <- as.factor(df3$individual)
  df4$individual <- as.factor(df4$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue",linetype = 1,size=0.8)+
    geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red",linetype = 1,size=0.8)+geom_line(aes_string(df4[,1],df4[,j],by=df4[,ncol(df4)]),colour="yellow",linetype = 1,size=0.8)+
      xlab("iteration") +scale_x_log10(breaks= c(10,100,200))+ ylab(expression(paste(E(V[pop]))))  
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}

c <- graphConvMC_se1(err_rwm_scaled[,c(1,3,9)],err_mix_scaled[,c(1,3,9)],err_mix_scaled[,c(1,3,9)],err_mix_scaled[,c(1,3,9)])
c <- graphConvMC_se1(err_rwm_scaled[,c(1,3,9)],err_mix_scaled[,c(1,3,9)],err_mix50_gamma2[,c(1,3,9)],err_mix50_gamma20[,c(1,3,9)])
d <- graphConvMC_se1(err_rwm_scaled[,c(1,6,9)],err_mix_scaled[,c(1,6,9)],err_mix50_gamma2[,c(1,6,9)],err_mix50_gamma20[,c(1,6,9)])


save <- grid.arrange(c,d, ncol=2)