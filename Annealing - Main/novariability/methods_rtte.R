setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Annealing - Main/novariability/R")
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
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Annealing - Main/novariability/")
source('mixtureFunctions.R') 
source('plots.R') 

library('rCMA')

###rtte
###RTTE
timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/Annealing - Main/novariability/data/rtte_data.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
saemix.data_rtte<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
timetoevent.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2]
  N <- nrow(psi)
  Nj <- length(T)
  censoringtime = 20
  lambda <- psi[id,1]
  beta <- psi[id,2]
  init <- which(T==0)
  cens <- which(T==censoringtime)
  ind <- setdiff(1:Nj, append(init,cens))
  hazard <- (beta/lambda)*(T/lambda)^(beta-1)
  H <- (T/lambda)^beta
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens] + H[cens-1]
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
  return(logpdf)
}

saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

saemix.model_rttenovar<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,0),ncol=2, 
  byrow=TRUE))

##RUNS

K1 = 150
K2 = 20
iterations = 1:(K1+K2+1)
end = K1+K2


#with var no sa
options_rtte<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),
  nbiter.saemix = c(K1,K2), nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0,an=FALSE,coeff=1, av=0)
rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
rtte <- cbind(iterations, rtte)


#RWM
replicate <- 3
final_optim <- 0
final_av <- 0
final_avnew <- 0

final_annealing <- 0
final_bayes <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(2,1),c(3,3),c(1,0.5),c(2.5,1.5))
  
 saemix.model_rttenovar<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(l[[m]],ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

  #No var
  ##### Optim (fmin search)
  options_rtte_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0),
    nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,an=FALSE,coeff=1, av=0)
  rtte_optim<-data.frame(saemix(saemix.model_rttenovar,saemix.data_rtte,options_rtte_with))
  rtte_optim <- cbind(iterations, rtte_optim)
  rtte_optim['individual'] <- m
  final_optim <- rbind(final_optim,rtte_optim)

  ##### AV
  options_rtte_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0),
    nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1,av=1)
  rtte_withav<-data.frame(saemix(saemix.model_rttenovar,saemix.data_rtte,options_rtte_with))
  rtte_withav <- cbind(iterations, rtte_withav)
  rtte_withav['individual'] <- m
  final_av <- rbind(final_av,rtte_withav)

  # ##### AV and new kernel
  # options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(1:3), av=1)
  # rtte_newkernelav<-data.frame(saemix(saemix.model_rttenovar,saemix.data_rtte,options_newkernel))
  # rtte_newkernelav <- cbind(iterations, rtte_newkernelav)
  # rtte_newkernelav['individual'] <- m
  # final_avnew <- rbind(final_avnew,rtte_newkernelav)

   ##### AV and new kernel
  options_an<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0),
    nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(0),an=TRUE,coeff=0.05, av=1)
  rtte_an<-data.frame(saemix(saemix.model_rttenovar,saemix.data_rtte,options_an))
  rtte_an <- cbind(iterations, rtte_an)
  rtte_an['individual'] <- m
  final_annealing <- rbind(final_annealing,rtte_an)

  ##### pseudo bayesian
  options_rtte_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,2),
    nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,an=FALSE,coeff=1,
     av=0,map.range=c(0))
  rtte_bayes<-data.frame(saemix(saemix.model_rttenovar,saemix.data_rtte,options_rtte_with))
  rtte_bayes <- cbind(iterations, rtte_bayes)
  rtte_bayes['individual'] <- m
  final_bayes <- rbind(final_bayes,rtte_bayes)
}

graphConvMC_diff4(final_optim,final_av,final_annealing,final_bayes, title="")

#black: optim
#blue: av
#red: annealing
#green: bayes

graphConvMC_diff(final_optim,final_annealing, title="")
graphConvMC_diff(final_optim,final_bayes, title="")