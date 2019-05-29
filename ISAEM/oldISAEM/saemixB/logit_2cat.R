
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/saemixB/R")
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

setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/saemixB/")
source("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/saemixB/plots.R")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################

# catModel <- inlineModel("
# [LONGITUDINAL]
# input =  {beta0,gamma0,delta0, dose}
# dose = {use=regressor}
# EQUATION:
# lm0 = beta0+gamma0*t + delta0*dose

# D = exp(lm0)+1

# p0 = exp(lm0)/D
# p1 = 1/D

# DEFINITION:
# y = {type=categorical, categories={0, 1}, 
#      P(y=0)=p0,
#      P(y=1)=p1}

# [INDIVIDUAL]
# input={beta0_pop, o_beta0,
#       gamma0_pop, o_gamma0,
#       delta0_pop, o_delta0}


# DEFINITION:
# beta0  ={distribution=normal, prediction=beta0_pop,  sd=o_beta0}

# gamma0  ={distribution=normal, prediction=gamma0_pop,  sd=o_gamma0}

# delta0  ={distribution=normal, prediction=delta0_pop,  sd=o_delta0}

# ")


# nobs = 10
# t<- seq(1, 100, by=nobs)
# reg <- list(name='dose',
#             time=t,
#             value=20)

# out  <- list(name='y', time=seq(1, 100, by=nobs))
# N  <- 200
# p <- c(beta0_pop=1, o_beta0=0, 
#        gamma0_pop= -1, o_gamma0=0.5,
#        delta0_pop=3, o_delta0=0.4)

# g1 <- list(size=N, parameter=p)
# res <- simulx(model=catModel, regressor = reg, output=out, group=g1)
# plot1 <- catplotmlx(res$y)
# print(plot1)

catModel <- inlineModel("
[LONGITUDINAL]
input =  {beta0,gamma0,delta0, dose}
dose = {use=regressor}
EQUATION:
lm0 = beta0+gamma0*t + delta0*dose

D = exp(lm0)+1
p0 = exp(lm0)/D
p1 = 1/D

DEFINITION:
y = {type=categorical, categories={0, 1}, 
     P(y=0)=p0,
     P(y=1)=p1}

[INDIVIDUAL]
input={beta0_pop, o_beta0,
      gamma0_pop, o_gamma0,
      delta0_pop, o_delta0}


DEFINITION:
beta0  ={distribution=normal, prediction=beta0_pop,  sd=o_beta0}
gamma0  ={distribution=normal, prediction=gamma0_pop,  sd=o_gamma0}
delta0  ={distribution=normal, prediction=delta0_pop,  sd=o_delta0}
")


nobs = 10
tobs<- seq(-20, 50, by=nobs)

reg1 <- list(name='dose',
            time=tobs,
            value=10*(tobs>0))

reg2 <- list(name='dose',
            time=tobs,
            value=20*(tobs>0))

reg3 <- list(name='dose',
            time=tobs,
            value=30*(tobs>0))

out  <- list(name='y', time=tobs)
N  <- 100
p <- c(beta0_pop= -4, o_beta0=0.3, 
       gamma0_pop= -1, o_gamma0=0.4,
       delta0_pop=1, o_delta0=0.3)

g1 <- list(size=N,regressor = reg1)
g2 <- list(size=N,regressor = reg2)
g3 <- list(size=N,regressor = reg3)
g <- list(g1,g2,g3)
res <- simulx(model=catModel,output=out, group=g,parameter=p)
plot1 <- catplotmlx(res$y)
print(plot1)


writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/saemixB/data/cat_data_isaem.csv")
cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/ISAEM/saemixB/data/cat_data_isaem.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"), name.predictors=c("y","dose","time"))



cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]
dose<-xidep[,2]
time<-xidep[,3]

beta0 <- psi[id,1]

gamma0 <- psi[id,2]

delta0 <- psi[id,3]

lm0 <- beta0+gamma0*time + delta0*dose

D <- exp(lm0)+1

P0 <- exp(lm0)/D
P1 <- 1/D

P.obs = log((level==0)*P0+(level==1)*P1)

return(P.obs)
}


cov <- matrix(c(1,0,0,
                0,1,0,
                0,0,1),ncol=3, byrow=TRUE)

saemix.model<-saemixModel(model=cat_data.model,description="cat model",type="likelihood",   
  psi0=matrix(c(-3,1,3),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("beta0",
    "gamma0",
    "delta0"))), 
  transform.par=c(0,0,0), fixed.estim=c(1,1,1),covariance.model=cov,omega.init=cov,error.model="constant")



K1 = 500
K2 = 100

iterations = 0:(K1+K2-1)
end = K1+K2
seed0 = 39546
map_range = c(1:6)

options<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,2),map.range=map_range,
  nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, nb.replacement=100,sampling='seq')
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])



options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,2),map.range=map_range, 
                          nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,2),map.range=map_range, nbiter.saemix = c(K1,K2),displayProgress=FALSE,
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])



# options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
#   nbiter.mcmc = c(2,2,2,2),map.range=map_range, nbiter.saemix = c(K1,K2),displayProgress=FALSE,
#   nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
# theo_mix75<-saemix_incremental(saemix.model,saemix.data,options.incremental75)
# theo_mix75 <- data.frame(theo_mix75$param)
# theo_mix75 <- cbind(iterations, theo_mix75[-1,])


options.incremental85<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,2),map.range=map_range, nbiter.saemix = c(K1,K2),displayProgress=FALSE,
  nbiter.sa=0,nbiter.burn =0, nb.replacement=85,sampling='randompass')
theo_mix85<-saemix_incremental(saemix.model,saemix.data,options.incremental85)
theo_mix85 <- data.frame(theo_mix85$param)
theo_mix85 <- cbind(iterations, theo_mix85[-1,])


theo_ref_scaled <- theo_ref
theo_mix50_scaled <- theo_mix50
theo_mix25_scaled <- theo_mix25
# theo_mix75_scaled <- theo_mix75
theo_mix85_scaled <- theo_mix85


theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25
# theo_mix75_scaled$iterations = theo_mix75_scaled$iterations*0.75
theo_mix85_scaled$iterations = theo_mix85_scaled$iterations*0.85

# graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix50_scaled,theo_mix50_scaled)
graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix85_scaled,theo_mix85_scaled)
# graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix85_scaled,theo_mix85_scaled)

