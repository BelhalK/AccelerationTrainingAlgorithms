####20 CHAINS
# load("warfa_isaem.RData")

# load("Rdata/pk100_1chain.RData")
# load("Rdata/pk100_10chains.RData")
# load("Rdata/pk100_20chains.RData")

source('R/aaa_generics.R') 
source('R/compute_LL.R') 
source('R/func_aux.R') 
source('R/func_distcond.R') 
source('R/func_FIM.R')
source('R/func_plots.R') 
source('R/func_simulations.R') 

source('R/main.R')
source('R/main_estep.R')
source('R/main_initialiseMainAlgo.R') 
source('R/main_mstep.R') 
source('R/SaemixData.R')
source('R/SaemixModel.R') 
source('R/SaemixRes.R') 

source('R/SaemixObject.R') 
source('R/zzz.R') 

source('R/main_incremental.R')
source('R/main_estep_incremental.R')

source('/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/incremental/R/mixtureFunctions.R')
source("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/incremental/plots.R")


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

N=500

param   <- c(
  ka_pop  = 2,    omega_ka  = 0.3,
  V_pop   = 10,   omega_V   = 0.2,
  Cl_pop  = 1,    omega_Cl  = 0.3, a =1)
  
res <- simulx(model     = model,
              parameter = param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(1,5,by=1)))

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

nchains = 50
gamma = 1
options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains,nbiter.mcmc = c(2,2,2,0), 
  nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, 
  map.range=c(0), nb.replacement=100,sampling='randomiter',gamma=gamma, algo="full")
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])

# graphConvMC_threekernels(theo_ref,theo_ref,theo_ref)
options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randomiter',gamma=gamma,algo="minibatch")
theo_mix75<-saemix_incremental(saemix.model,saemix.data,options.incremental75)
theo_mix75 <- data.frame(theo_mix75$param)
theo_mix75 <- cbind(iterations, theo_mix75[-1,])


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='seq',gamma=gamma,algo="minibatch")
theo_mix50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo_mix50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = nchains, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='seq',gamma=gamma,algo="minibatch")
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])

# graphConvMC_threekernels(theo_ref,theo_mix25,theo_mix50)

theo_ref_scaled <- theo_ref
theo_mix25_scaled <- theo_mix25
theo_mix50_scaled <- theo_mix50
theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
# graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled)
theo_mix75_scaled <- theo_mix75
theo_mix75_scaled$iterations = theo_mix75_scaled$iterations*0.75

graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix50_scaled,theo_mix75_scaled)



theo_ref_scaled$algo <- 'full'
theo_mix25_scaled$algo <- 'quarter'
theo_mix50_scaled$algo <- 'half'
theo_mix75_scaled $algo<- 'three quarter'
comparison <- rbind(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix75_scaled)
var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

# write.csv(var, file = "../notebooks/data/pk100_1chain.csv")
# write.csv(var, file = "../notebooks/data/pk100_5chains.csv")
# write.csv(var, file = "../notebooks/data/pk100_10chains.csv")
# write.csv(var, file = "../notebooks/data/pk100_20chains.csv")

# save.image("Rdata/pk100_1chain.RData")
# save.image("Rdata/pk100_5chains.RData")
# save.image("Rdata/pk100_10chains.RData")
# save.image("Rdata/pk100_20chains.RData")


graphConvMC_threekernels <- function(df,df2,df3, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red")+
      xlab("iteration") + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}
