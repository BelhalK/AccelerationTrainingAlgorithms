# save.image("zifro_novar_an.RData")
load("zifro_novar_an.RData")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/Annealing - Main/novariability/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('main.R')
  source('main_estep.R')
  source('main_initialiseMainAlgoold.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/Annealing - Main/novariability/")
source('plots.R') 
source('mixtureFunctions.R') 

library('rCMA')
library(mlxR)
###zifro


# zifro_data <- read.csv("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/Annealing - Main/novariability/data/zifro.csv", header=T,sep=",")
# saemix.data_zifro<-saemixData(name.data=zifro_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("amount","time"),name.response=c("y"), name.X="x")

zifro_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/Annealing - Main/novariability/data/dataPK_zifrosilone.txt", header=T)
saemix.data_zifro<-saemixData(name.data=zifro_data,header=TRUE,sep=" ",na=NA, name.group=c("ID"),
  name.predictors=c("AMT","TIME"),name.response=c("Y"), name.X="X")



# model <- inlineModel("
#               [INDIVIDUAL]
#               input = {Tlag_pop, omega_Tlag, ka_pop, omega_ka, V_pop, omega_V, alpha_pop ,omega_alpha, beta_pop, omega_beta}

#               DEFINITION:
#               Tlag = {distribution=lognormal, prediction=Tlag_pop, sd=omega_Tlag}
#               ka = {distribution=lognormal, prediction=ka_pop, sd=omega_ka}
#               V = {distribution=lognormal, prediction=V_pop,sd=omega_V}
#               alpha = {distribution=lognormal, prediction=alpha_pop,sd=omega_alpha}
#               beta = {distribution=lognormal, prediction=beta_pop,sd=omega_beta}


#               [LONGITUDINAL]
#               input = {Tlag, ka, V, alpha, beta,a}

#               EQUATION:
#               Cc = pkmodel(Tlag, ka, V, Cl=alpha*(V^beta))

#               OUTPUT:
#               output = Cc

#               DEFINITION:
#               y1 = {distribution=normal, prediction=Cc, errorModel=constant(a)}

#                       ")

# adm  <- list(amount=100, time=seq(0,50,by=50))


# p <- c(Tlag_pop=0.000359, omega_Tlag=3.5,  
#       ka_pop=0.7234, omega_ka=0.8,
#        V_pop=195, omega_V=0.6, 
#        alpha_pop=1.3, omega_alpha=0,  
#        beta_pop=0.6, omega_beta=0,  
#        a=0.2)
# y1 <- list(name='y1', time=seq(1,to=50,by=5))

# res <- simulx(model = model,
#                  treatment = adm,
#                  parameter = p,
#                  group = list(size=15, level="individual"),
#                  output = y1)

# zifro_data <- res$y1
# zifro_data$amount = 100
# head(zifro_data)
# saemix.data_zifro<-saemixData(name.data=zifro_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("amount","time"),name.response=c("y1"), name.X="x")

saemix.data_zifro<-saemixData(name.data=zifro_data,header=TRUE,sep=" ",na=NA, name.group=c("ID"),
  name.predictors=c("AMT","TIME"),name.response=c("Y"), name.X="x")


model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
  Tlag<-psi[id,1]
  ka<-psi[id,2]
  V<-psi[id,3]
  alpha<-psi[id,4]
  beta<-psi[id,5]
  CL<-alpha*V^beta
  k<-CL/V
  # ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  dt <- pmax(time-Tlag, 0)
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*dt)-exp(-ka*dt))
  return(ypred)
}


saemix.model_zifro<-saemixModel(model=model1cpt,description="zifrorin",type="structural"
  ,psi0=matrix(c(0.2,1,500,1,1),ncol=5,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","alpha","beta"))),
  transform.par=c(1,1,1,1,1),omega.init=matrix(c(0.1,0,0,0,0,0,0.1,0,0,0,0,0,0.1,0,0,0,0,0,0.1,0,0,0,0,0,0.1),ncol=5,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5, 
  byrow=TRUE),error.model="constant")

# saemix.model_zifronovar<-saemixModel(model=model1cpt,description="zifrorin"
#   ,psi0=matrix(c(0.2,1,250,1,1),ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
#   transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol=5, 
#   byrow=TRUE),error.model="constant")

saemix.model_zifronovar<-saemixModel(model=model1cpt,description="zifrorin",type="structural"
  ,psi0=matrix(c(0.158,2,200,0.5,2),ncol=5,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","alpha","beta"))),
  transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol=5, 
  byrow=TRUE),fixed.estim=c(1,1,1,1,1),error.model="constant")




K1 = 300
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2


#With var 
options_zifro_without<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
zifro_without<-data.frame(saemix(saemix.model_zifro,saemix.data_zifro,options_zifro_without))
zifro_without <- cbind(iterations, zifro_without)

#no var and optim
options_zifro_without<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0,0), nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
zifro_without<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_without))
zifro_without <- cbind(iterations, zifro_without)

graphConvMC_twokernels(zifro_without,zifro_without)

#no var but av 
options_zifro_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0), nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=1)
zifro_without<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_without))
zifro_without <- cbind(iterations, zifro_without)


#No var but pseudobayes
options_zifro_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
zifro_withnosa<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_with))
zifro_withnosa <- cbind(iterations, zifro_withnosa)
graphConvMC_twokernels(zifro_without[,1:6],zifro_withnosa)


replicate = 3
seed0 = 395246

#RWM
final_optim <- 0
final_av <- 0
final_avnew <- 0
final_bayes <- 0

final_annealing <- 0


for (m in 1:replicate){
  print(m)
  l = list(c(0.1,0.18,100,1,1),c(0.15,0.2,150,1.5,1.5),c(0.11,0.15,130,1.2,1.2),c(0.2,0.18,210,1,1),c(0.25,0.18,50,1,1))
  
# saemix.model_zifro<-saemixModel(model=model1cpt,description="zifrorin",type="structural"
#   ,psi0=matrix(l[[m]],ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
#   transform.par=c(1,1,1,1,1),omega.init=matrix(c(0.1,0,0,0,0,0,0.1,0,0,0,0,0,0.1,0,0,0,0,0,0.1,0,0,0,0,0,0.1),ncol=5,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5, 
#   byrow=TRUE),error.model="constant")

saemix.model_zifronovar<-saemixModel(model=model1cpt,description="zifrorin",type="structural"
  ,psi0=matrix(l[[m]],ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
  transform.par=c(1,1,1,1,1),
  omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol=5, 
  byrow=TRUE),error.model="constant")


  # #No var
  # ##### Optim (fmin search)
  # options_zifro_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0),
  #   nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,an=FALSE,coeff=1,nbiter.burn =0, av=0)
  # zifro_optim<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_with))
  # zifro_optim <- cbind(iterations, zifro_optim)
  # zifro_optim['individual'] <- m
  # final_optim <- rbind(final_optim,zifro_optim)
  
  # #### AV
  # options_zifro_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0),
  #   nbiter.sa=K1,nbiter.saemix = c(K1,K2),displayProgress=FALSE,an=FALSE,coeff=1,nbiter.burn =0, map.range=c(0),av=1)
  # zifro_withav<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_with))
  # zifro_withav <- cbind(iterations, zifro_withav)
  # zifro_withav['individual'] <- m
  # final_av <- rbind(final_av,zifro_withav)

  # ##### AV and new kernel
  # options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6,0),
  # nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(1), av=1)
  # zifro_newkernelav<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_newkernel))
  # zifro_newkernelav <- cbind(iterations, zifro_newkernelav)
  # zifro_newkernelav['individual'] <- m
  # final_avnew <- rbind(final_avnew,zifro_newkernelav)


  # ##### Annealing MCMC
  options_an<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0),
  nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, an=TRUE,coeff=0.5,
  map.range=c(0), av=1)
  zifro_an<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_an))
  zifro_an <- cbind(iterations, zifro_an)
  zifro_an['individual'] <- m
  final_annealing<- rbind(final_annealing,zifro_an)

  ##### pseudo bayesian

  # options_zifro_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,2),
  #   nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,an=FALSE,coeff=1, 
  #   nbiter.burn =0, av=0,map.range=c(0))
  # zifro_bayes<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_with))
  # zifro_bayes <- cbind(iterations, zifro_bayes)
  # zifro_bayes['individual'] <- m
  # final_bayes <- rbind(final_bayes,zifro_bayes)

}


graphConvMC_diff4(final_optim,final_av,final_annealing,final_bayes, title="")
#black: optim
#blue: av
#red: annealing
#green: bayes

# graphConvMC_diff4(final_optim,final_av,final_avnew,final_bayes, title="")
