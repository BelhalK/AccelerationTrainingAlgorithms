setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/novariability/R")
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
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/novariability/")
source('plots.R') 
source('mixtureFunctions.R') 

library('rCMA')
###zifro


# zifro_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/data/dataPK_zifrosilone.csv", header=T,sep=";")
zifro_data <- read.table("/Users/karimimohammedbelhal/Desktop/novariability/data/dataPK_zifrosilone.txt", header=T)
saemix.data_zifro<-saemixData(name.data=zifro_data,header=TRUE,sep=" ",na=NA, name.group=c("ID"),
  name.predictors=c("AMT","TIME"),name.response=c("Y"), name.X="X")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  T<-psi[id,1]
  ka<-psi[id,2]
  V<-psi[id,3]
  alpha<-psi[id,4]
  beta<-psi[id,5]
  CL<-alpha*V^beta
  k<-CL/V
  # ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*(tim-T))-exp(-ka*(tim-T)))
  return(ypred)
}


saemix.model_zifro<-saemixModel(model=model1cpt,description="zifrorin",type="structural"
  ,psi0=matrix(c(0.2,1,250,1,1),ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
  transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5, 
  byrow=TRUE),error.model="exponential")

# saemix.model_zifronovar<-saemixModel(model=model1cpt,description="zifrorin"
#   ,psi0=matrix(c(0.2,1,250,1,1),ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
#   transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol=5, 
#   byrow=TRUE),error.model="exponential")

saemix.model_zifronovar<-saemixModel(model=model1cpt,description="zifrorin",type="structural"
  ,psi0=matrix(c(0.158,0.18,40,1,1),ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
  transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol=5, 
  byrow=TRUE),error.model="exponential")


K1 = 400
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2



#With var 
options_zifro_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
zifro_without<-data.frame(saemix(saemix.model_zifro,saemix.data_zifro,options_zifro_without))
zifro_without <- cbind(iterations, zifro_without)

#no var
options_zifro_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0), nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
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
for (m in 1:replicate){
  print(m)
  l = list(c(0.158,0.18,40,1,1),c(0.168,0.2,43,1,1),c(0.138,0.15,37,1,1),c(0.158,0.18,40,1,1))
  

saemix.model_zifro<-saemixModel(model=model1cpt,description="zifrorin",type="structural"
  ,psi0=matrix(l[[m]],ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
  transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5, 
  byrow=TRUE),error.model="exponential")

saemix.model_zifronovar<-saemixModel(model=model1cpt,description="zifrorin",type="structural"
  ,psi0=matrix(l[[m]],ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
  transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol=5, 
  byrow=TRUE),error.model="exponential")


  #No var
  ##### Optim (fmin search)
  options_zifro_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
  zifro_optim<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_with))
  zifro_optim <- cbind(iterations, zifro_optim)
  zifro_optim['individual'] <- m
  final_optim <- rbind(final_optim,zifro_optim)
  ##### AV
  options_zifro_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0),nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, map.range=c(0),av=1)
  zifro_withav<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_with))
  zifro_withav <- cbind(iterations, zifro_withav)
  zifro_withav['individual'] <- m
  final_av <- rbind(final_av,zifro_withav)

  ##### AV and new kernel
  options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6,0),nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1), av=1)
  zifro_newkernelav<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_newkernel))
  zifro_newkernelav <- cbind(iterations, zifro_newkernelav)
  zifro_newkernelav['individual'] <- m
  final_avnew <- rbind(final_avnew,zifro_newkernelav)

  ##### pseudo bayesian

  options_zifro_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0,map.range=c(0))
  zifro_bayes<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_with))
  zifro_bayes <- cbind(iterations, zifro_bayes)
  zifro_bayes['individual'] <- m
  final_bayes <- rbind(final_bayes,zifro_bayes)
}



graphConvMC_diff4(final_optim,final_av,final_avnew,final_bayes, title="")

#black: optim
#blue: av
#red: av newkernel
#green: bayes


