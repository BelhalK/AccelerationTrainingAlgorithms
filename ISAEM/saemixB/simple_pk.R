
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
  time<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  k <- Cl/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*time)-exp(-ka*time))
  return(ypred)
}


ka_true <- 2
V_true <- 13
Cl_true <- 1
o_ka <- sqrt(0.2)
o_V <- sqrt(0.1)
o_Cl <- sqrt(0.3)
a_true<-1



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

N=20000
pop.param   <- c(
  ka_pop  = ka_true,    omega_ka  = o_ka,
  V_pop   = V_true,   omega_V   = o_V,
  Cl_pop  = Cl_true,    omega_Cl  = o_Cl, a =a_true)
  
res <- simulx(model     = myModel,
              parameter = pop.param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(0,4,by=1)))
  
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



K1 = 300
K2 = 100
iterations = 0:(K1+K2-1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50

seed0=3456


options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='seq')
theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref[-1,])


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo_mix50<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental50))
theo_mix50 <- cbind(iterations, theo_mix50[-1,])



options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_mix25[-1,])

graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix25_scaled,theo_mix25_scaled)

options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
theo_mix75<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental75))
theo_mix75 <- cbind(iterations, theo_mix75[-1,])


options.incremental85<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
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

graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix75_scaled,theo_mix85_scaled)
