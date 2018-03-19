# savedwarfarin.saemix <- warfarin.saemix
# save.image("pkcov.RData")
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


model2 <- inlineModel("
                      [LONGITUDINAL]
                      input = {base, slope,a}
                      EQUATION:
                      Cc = base+slope*t

                      DEFINITION:
                      y1 ={distribution=normal, prediction=Cc, sd=a}
                      
                      [INDIVIDUAL]
                      input={base_pop,o_base, slope_pop,o_slope}
                      
                      DEFINITION:
                      base  ={distribution=lognormal, prediction=base_pop,  sd=o_base}
                      slope  ={distribution=normal, prediction=slope_pop,  sd=o_slope}
                      ")


p <- c(base_pop=140, o_base=0.3,
       slope_pop=1, o_slope=0.4,  
       a=1)
y1 <- list(name='y1', time=seq(1,to=50,by=5))


res2a2 <- simulx(model = model2,
                 parameter = p,
                 group = list(size=100, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/oxboys/ox_synth.csv")

oxboys.saemix<-read.table( "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/oxboys/ox_synth.csv",header=T,na=".",sep=",")
oxboys.saemix_less <- oxboys.saemix[,]
n <- length(unique(oxboys.saemix_less$id))

saemix.data<-saemixData(name.data=oxboys.saemix_less,header=TRUE,
  name.group=c("id"),name.predictors=c("time"),name.response=c("y"),
  units=list(x="yr",y="cm"))


growth.linear<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (2 columns, base and slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  base<-psi[id,1]
  slope<-psi[id,2]
  f<-base+slope*x
  return(f)
}

saemix.model<-saemixModel(model=growth.linear,description="Linear model",type="structural",
  psi0=matrix(c(110,5),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(1,0),covariance.model=matrix(c(1,0,0,1),ncol=2,byrow=TRUE),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE), 
  error.model="constant")


##RUNS

K1 = 100
K2 = 50
iterations = 1:(K1+K2+1)
end = K1+K2

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
warfa<-data.frame(saemix(saemix.model,saemix.data,options_warfa))
warfa <- cbind(iterations, warfa)


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(1:3))
warfanew<-data.frame(saemix(saemix.model,saemix.data,options_warfanew))
warfanew <- cbind(iterations, warfanew)


graphConvMC_twokernels(warfa,warfanew)

options.new<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,6), nbiter.saemix = c(K1,K2)
  ,nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(1:3), nb.replacement=10,sampling='randompass')
theo_new<-data.frame(saemix_incremental(saemix.model,saemix.data,options.new))
theo_new <- cbind(iterations, theo_new)

graphConvMC_twokernels(warfanew,theo_new)

#isaem rwm
options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='seq')
theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo_mix50<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental50))
theo_mix50 <- cbind(iterations, theo_mix50)



options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_mix25)


options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
theo_mix75<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental75))
theo_mix75 <- cbind(iterations, theo_mix75)


options.incremental85<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=85,sampling='randompass')
theo_mix85<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental85))
theo_mix85 <- cbind(iterations, theo_mix85)


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


#isaem new kernel
options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,6), nbiter.saemix = c(K1,K2)
  ,nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(1:3), nb.replacement=100,sampling='seq')
theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,6), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:3),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo_mix50<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental50))
theo_mix50 <- cbind(iterations, theo_mix50)



options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,6), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:3),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_mix25)


options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,6), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:3),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
theo_mix75<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental75))
theo_mix75 <- cbind(iterations, theo_mix75)


options.incremental85<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,6), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:3),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=85,sampling='randompass')
theo_mix85<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental85))
theo_mix85 <- cbind(iterations, theo_mix85)


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




