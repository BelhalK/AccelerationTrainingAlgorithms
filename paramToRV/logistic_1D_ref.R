library(saemix)

logit_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/data/Data_logistic_1D_simule.txt", header=T)
saemix.data_logit<-saemixData(name.data=logit_data,header=TRUE,sep=" ",na=NA, name.group=c("group"),
  name.predictors=c("X"),name.response=c("Y"), name.X="X")

model1cpt<-function(psi,id,xidep) { 
  tim<-xidep[,1]  
  p0<-psi[id,1]
  alpha<-psi[id,2]
  tau<-psi[id,3]
  ypred<- 1/(1+((1/p0)-1)*exp(-alpha*(tim-tau)/(p0*(1-p0))))
  return(ypred)
}



saemix.model_logit<-saemixModel(model=model1cpt,description="logitrin"
  ,psi0=matrix(c(0.6,0.05,60),ncol=3,byrow=TRUE, dimnames=list(NULL, c("p0","alpha","tau"))),
  transform.par=c(3,2,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

saemix.model_logitnovar<-saemixModel(model=model1cpt,description="logitrin"
  ,psi0=matrix(c(0.6,0.05,60),ncol=3,byrow=TRUE, dimnames=list(NULL, c("p0","alpha","tau"))),
  transform.par=c(3,2,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE))


K1 = 400
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2


#With var
options_logit_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0)
logit_without<-data.frame(saemix(saemix.model_logit,saemix.data_logit,options_logit_without))
logit_without <- cbind(iterations, logit_without)

#no var
options_logit_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0)
logit_without<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_logit_without))
logit_without <- cbind(iterations, logit_without)