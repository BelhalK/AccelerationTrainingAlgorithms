
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
  source('mixtureFunctions.R')

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/")

###logit


logit_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/Data_logistic_1D_simule.txt", header=T)
saemix.data_logit<-saemixData(name.data=logit_data,header=TRUE,sep=" ",na=NA, name.group=c("group"),
  name.predictors=c("X"),name.response=c("Y"), name.X="X")

model1cpt<-function(psi,id,xidep) { 
  tim<-xidep[,1]  
  p0<-psi[id,1]
  alpha<-psi[id,2]
  tau<-psi[id,3]
  ypred<-1/(1+((1/p0)-1)*exp(-alpha*(tim-tau)/(p0*(1-p0))))
  return(ypred)
}


saemix.model_logistic<-saemixModel(model=model1cpt,description="logitrin",type="structural"
  ,psi0=matrix(c(0.6,0.05,60),ncol=3,byrow=TRUE, dimnames=list(NULL, c("p0","alpha","tau"))),
  transform.par=c(3,1,0),omega.init=matrix(c(1,0,0,0,1,0,0,0,100),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))



K1 = 200
K2 = 10
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize<-50

options_logistic<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
logistic<-data.frame(saemix(saemix.model_logistic,saemix.data_logit,options_logistic))
logistic<-cbind(iterations,logistic)

options_logisticincr<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize)
logisticincr<-data.frame(saemix_incremental(saemix.model_logistic,saemix.data_logit,options_logisticincr))
logisticincr<-cbind(iterations,logisticincr)

graphConvMC2_saem(logistic,logisticnew, title="new kernel")


logistic$algo <- 'rwm'
logisticincr$algo <- 'ISAEM'

logistic_scaled <- logistic[rep(seq_len(nrow(logistic)), each=100/batchsize),]
logistic_scaled$iterations = 1:(2*(K1+K2+1))


comparison <- 0
# comparison <- rbind(theo_ref,theo_incremental)
comparison <- rbind(logistic_scaled[iterations,],logisticincr)

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)

options_logisticnew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5))
logisticnew<-data.frame(saemix(saemix.model_logistic,saemix.data_logit,options_logisticnew))
logisticnew<-cbind(iterations,logisticnew)
