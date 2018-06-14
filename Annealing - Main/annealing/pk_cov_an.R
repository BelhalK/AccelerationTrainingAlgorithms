library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
# save.image("pk_cov_an.RData")
load("pk_cov_an.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/annealing/R")
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

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/annealing")
source('plots.R') 


Tlag_true=0.78
ka_true <- 1
V_true <- 8
Cl_true <- 0.1

o_Tlag <- 0.57 #o^2=0.32
o_ka <- 0.5 #o^2=0.25
o_V <- 0.2  #o^2=0.04
o_Cl <- 0.3  #o^2=0.09
a_true = 0.266
beta_V_lw70_true = 0.8818
beta_Cl_lw70_true = 0.60411

  myModel <- inlineModel("

[COVARIATE]
input = wt

EQUATION:
lw70 = log(wt/70)

[INDIVIDUAL]
input = {Tlag_pop, omega_Tlag, ka_pop, omega_ka, V_pop, beta_V_lw70, lw70, omega_V, Cl_pop, beta_Cl_lw70, omega_Cl}

DEFINITION:
Tlag = {distribution=lognormal, typical=Tlag_pop, sd=omega_Tlag}
ka = {distribution=lognormal, typical=ka_pop, sd=omega_ka}
V = {distribution=lognormal, typical=V_pop, covariate=lw70, coefficient=beta_V_lw70, sd=omega_V}
Cl = {distribution=lognormal, typical=Cl_pop, covariate=lw70, coefficient=beta_Cl_lw70, sd=omega_Cl}

[LONGITUDINAL]
input =  {Tlag, ka, V, Cl,a}

EQUATION:
Cc = pkmodel(Tlag, ka, V, Cl)

OUTPUT:
output = {Cc}

DEFINITION:
y1 = {distribution=normal, prediction=Cc, sd=a}
")


N <- 30
populationParameter   <- c(Tlag_pop= Tlag_true, omega_Tlag= o_Tlag,
  ka_pop  = ka_true,    omega_ka  = o_ka,
  V_pop   = V_true,   omega_V   = o_V,
  Cl_pop  = Cl_true,    omega_Cl  = o_Cl, a =a_true, beta_V_lw70 = beta_V_lw70_true, beta_Cl_lw70 = beta_Cl_lw70_true)

individualCovariate <- data.frame(matrix(NA, nrow = N, ncol = 2))
colnames(individualCovariate) <- c("id","wt")
individualCovariate$id <- 1:N
individualCovariate$wt <- runif(N, 60.0, 80.0)
list.param <- list(populationParameter,individualCovariate)

amount <- 100
res <- simulx(model     = myModel,
              parameter = list.param,
              treatment = list(time=0, amount=amount),
              output    = list(name='y1', time=seq(1,3,by=1)))


# call the simulator 
warfarin.saemix <- res$y1
individualCovariate$wt <- log(individualCovariate$wt/70)
warfarin.saemix$amount <- amount
warfarin.saemix <- merge(individualCovariate ,warfarin.saemix,by="id")

saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time", name.covariates=c("wt"))

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
  Tlag<-psi[id,1]
  ka<-psi[id,2]
  V<-psi[id,3]
  Cl<-psi[id,4]
  k<-Cl/V
  dt <- pmax(time-Tlag, 0)
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*dt)-exp(-ka*dt))
  return(ypred)
}

saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.2,3,10,2),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","Cl"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE),covariate.model=matrix(c(0,0,1,1),ncol=4,byrow=TRUE),error.model="constant")


K1 = 300
K2 = 100
iterations = 1:(K1+K2)
end = K1+K2

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1)
warfa<-data.frame(saemix(saemix.model,saemix.data,options_warfa))
warfa <- cbind(iterations, warfa[-1,])

options_warfa.sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1)
warfa.sa<-data.frame(saemix(saemix.model,saemix.data,options_warfa.sa))
warfa.sa <- cbind(iterations, warfa.sa[-1,])


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, an=TRUE,coeff=0.0303)
warfanew<-data.frame(saemix(saemix.model,saemix.data,options_warfanew))
warfanew <- cbind(iterations, warfanew[-1,])

plot3(warfa,warfa.sa,warfanew)

# options_warfanew2<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, an=TRUE,coeff=1)
# warfanew2<-data.frame(saemix(saemix.model,saemix.data,options_warfanew2))
# warfanew2 <- cbind(iterations, warfanew2[-1,])

# options_warfanew3<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, an=TRUE,coeff=2)
# warfanew3<-data.frame(saemix(saemix.model,saemix.data,options_warfanew3))
# warfanew3 <- cbind(iterations, warfanew3[-1,])



plot5(warfa,warfa.sa,warfanew,warfanew2,warfanew3)



replicate = 3

final.ref <- 0
final.sa <- 0
final.an <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(0.2,3,10,2),c(0.5,5,12,1),c(0.7,6,14,4))
  
  saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(l[[m]],ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","Cl"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE),covariate.model=matrix(c(0,0,1,1),ncol=4,byrow=TRUE),error.model="constant")



  options<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1)
  pk<-data.frame(saemix(saemix.model,saemix.data,options))
  pk <- cbind(iterations, pk[-1,])
  pk['individual'] <- m
  final.ref <- rbind(final.ref,pk)


  options.sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1, alpha.sa=0.95)
  pk.sa<-data.frame(saemix(saemix.model,saemix.data,options.sa))
  pk.sa <- cbind(iterations, pk.sa[-1,])
  pk.sa['individual'] <- m
  final.sa <- rbind(final.sa,pk.sa)

  optionsnew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, an=TRUE,coeff=0.0303)
  pknew<-data.frame(saemix(saemix.model,saemix.data,optionsnew))
  pknew <- cbind(iterations, pknew[-1,])
  pknew['individual'] <- m
  final.an <- rbind(final.an,pknew)

}


diff(final.ref,final.sa,final.an)
