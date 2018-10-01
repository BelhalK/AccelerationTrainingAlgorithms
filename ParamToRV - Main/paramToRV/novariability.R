setwd("/Users/karimimohammedbelhal/Desktop/paramToRV/saemixnonvariability")
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
setwd("/Users/karimimohammedbelhal/Desktop/paramToRV/")
source('plots.R') 
###WARFA
warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/paramToRV/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]
  CL<-k*V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


###RTTE
timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/paramToRV/data/rtte_data.csv", header=T, sep=",")
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


##RUNS

K1 = 150
K2 = 20
iterations = 1:(K1+K2+1)
end = K1+K2

options_warfa_without<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_without<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_without))
warfa_without <- cbind(iterations, warfa_without)

graphConvMC_twokernels(warfa_without,warfa_without)

options_warfa_without_with_sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_without_sa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_without_with_sa))
warfa_without_sa <- cbind(iterations, warfa_without_sa)

graphConvMC_twokernels(warfa_without,warfa_without_sa)
graphConvMC_twokernelslog(warfa_without,warfa_without_sa)

options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=1)
warfa_withav<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_with))
warfa_withav <- cbind(iterations, warfa_withav)

graphConvMC_twokernels(warfa_withav,warfa_without)
graphConvMC_twokernelslog(warfa_withav,warfa_without_sa)

warfa_newkernel_sanovar <- warfa_newkernel_sa
graphConvMC_twokernelslog(warfa_newkernel_sa[,1:6],warfa_without_sa[,1:6])

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
warfa_newkernel<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_newkernel))
warfa_newkernel <- cbind(iterations, warfa_newkernel)
graphConvMC_twokernelslog(warfa_without,warfa_newkernel)
graphConvMC_twokernels(warfa_without,warfa_newkernel)

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
warfa_newkernel_sa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_newkernel))
warfa_newkernel_sa <- cbind(iterations, warfa_newkernel_sa)
graphConvMC_twokernelslog(warfa_without_sa,warfa_newkernel_sa)

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2), nbiter.sa=K1,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=1)
warfa_newkernelav<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_newkernel))
warfa_newkernelav <- cbind(iterations, warfa_newkernelav)

#AV=1
graphConvMC_twokernelslog(warfa_withav,warfa_newkernelav, title="AV")
#Simulate annealing
graphConvMC_twokernelslog(warfa_without_sa,warfa_newkernel_sa, title="SA")

#Weibull
options_rtte<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2), nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, av=0)
rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
rtte <- cbind(iterations, rtte)

options_rttesa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
rttesa<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttesa))
rttesa <- cbind(iterations, rttesa)

graphConvMC_twokernelslog(rtte,rttesa)

options_rtteav<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=1)
rtteav<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtteav))
rtteav <- cbind(iterations, rtteav)


options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.saemix = c(K1,K2), nbiter.sa=0,displayProgress=TRUE,map.range=c(1:5),nbiter.burn =0, av=0)
rttenew<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenew))
rttenew <- cbind(iterations, rttenew)

graphConvMC_twokernelslog(rtte,rttenew)

options_rttenewsa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.saemix = c(K1,K2),displayProgress=TRUE,map.range=c(1:5),nbiter.burn =0, av=0)
rttenewsa<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenewsa))
rttenewsa <- cbind(iterations, rttenewsa)

options_rttenewav<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.saemix = c(K1,K2),displayProgress=TRUE,map.range=c(1:5),nbiter.burn =0, av=1)
rttenewav<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenewav))
rttenewav <- cbind(iterations, rttenewav)

#AV=1
graphConvMC_twokernelslog(rtteav,rttenewav, title="AV")
#Simulate annealing
graphConvMC_twokernelslog(rttesa,rttenewsa, title="SA")



PD1.saemix<-read.table( "data/PD1.saemix.tab",header=T,na=".")
PD2.saemix<-read.table( "data/PD2.saemix.tab",header=T,na=".")
PD2.saemix <- subset(PD2.saemix, dose!="90")
saemix.data1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))
saemix.data2<-saemixData(name.data=PD2.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))

modelemax<-function(psi,id,xidep) {
# input:
# psi : matrix of parameters (3 columns, E0, Emax, EC50)
# id : vector of indices
# xidep : dependent variables (same nb of rows as length of id)
# returns:
# a vector of predictions of length equal to length of id
dose<-xidep[,1]
e0<-psi[id,1]
emax<-psi[id,2]
e50<-psi[id,3]
f<-e0+emax*dose/(e50+dose)
return(f)
}

saemix.model<-saemixModel(model=modelemax,description="Emax model",type="structural",
psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3,
byrow=TRUE),error.model="constant")

options_pd_without<-list(seed=39546,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, av=0)
pd_without<-data.frame(saemix(saemix.model,saemix.data1,options_pd_without))
pd_without <- cbind(iterations, pd_without)

options_pd_without_with_sa<-list(seed=39546,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
pd_without_with_sa<-data.frame(saemix(saemix.model,saemix.data1,options_pd_without_with_sa))
pd_without_with_sa <- cbind(iterations, pd_without_with_sa)

options_pd_with<-list(seed=39546,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=1)
pd_with<-data.frame(saemix(saemix.model,saemix.data1,options_pd_with))
pd_with <- cbind(iterations, pd_with)

# options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
# pd_newkernel<-data.frame(saemix(saemix.model,saemix.data2,options_newkernel))
# pd_newkernel <- cbind(iterations, pd_newkernel)
# graphConvMC_twokernels(pd_without,pd_newkernel)
# graphConvMC_twokernelslog(pd_without,pd_newkernel)

