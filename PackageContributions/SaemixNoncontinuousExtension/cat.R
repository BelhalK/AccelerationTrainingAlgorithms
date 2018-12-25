
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/PackageContributions/SaemixNoncontinuousExtension/R")
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
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/PackageContributions/SaemixNoncontinuousExtension")
library("mlxR")

cat_data.saemix<-read.table("data/logistic2cat_test.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"), name.predictors=c("y","dose","time"))


plot1 <- catplotmlx(res$y)
print(plot1)


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

P.obs = (level==0)*P0+(level==1)*P1

return(P.obs)
}


saemix.model<-saemixModel(model=cat_data.model,description="cat model",type="likelihood",   
  psi0=matrix(c(2,1,2),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))), 
  transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")


K1 = 500
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2


options<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)
cat.fit<-saemix(saemix.model,saemix.data,options)
