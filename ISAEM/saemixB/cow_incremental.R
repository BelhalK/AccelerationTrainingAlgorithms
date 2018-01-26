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
source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/plots.R")


require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

# Doc
# Doc
data(cow.saemix)
saemix.data<-saemixData(name.data=cow.saemix,header=TRUE,name.group=c("cow"), 
  name.predictors=c("time"),name.response=c("weight"), 
  name.covariates=c("birthyear","twin","birthrank"), 
  units=list(x="days",y="kg",covariates=c("yr","-","-")))


# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
growthcow<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, a, b, k)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  a<-psi[id,1]
  b<-psi[id,2]
  k<-psi[id,3]
  f<-a*(1-b*exp(-k*x))
  return(f)
}
saemix.model<-saemixModel(model=growthcow,type="structural",
  description="Exponential growth model", 
  psi0=matrix(c(700,0.9,0.02,0,0,0),ncol=3,byrow=TRUE, 
  dimnames=list(NULL,c("A","B","k"))),transform.par=c(1,1,1),fixed.estim=c(1,1,1), 
  covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
  omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")



K1 = 300
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50
replicate = 3
seed0 = 395246

#RWM
final_rwm <- 0
final_50 <- 0
final_25 <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(500,3,0.02),c(600,2,0.02),c(700,4,0.02),c(700,0.9,0.02))


saemix.model<-saemixModel(model=growthcow,type="structural",
  description="Exponential growth model", 
  psi0=matrix(l[[m]],ncol=3,byrow=TRUE, 
  dimnames=list(NULL,c("A","B","k"))),transform.par=c(1,1,1),fixed.estim=c(1,1,1), 
  covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
  omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE),error.model="constant")

  options_cow<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
  cow<-data.frame(saemix(saemix.model,saemix.data,options_cow))
  cow<-cbind(iterations,cow)
  cow['individual'] <- m
  cow$algo <- 'rwm'
  cow_scaled <- cow[-1,]
  cow_scaled$iterations = seq(1, 4*end, by=4)
  cow_scaled <- cow_scaled[rep(seq_len(nrow(cow_scaled)), each=4),]
  final_rwm <- rbind(final_rwm,cow_scaled[0:end,])

  options_cowincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
  cowincr25<-data.frame(saemix_incremental(saemix.model,saemix.data,options_cowincr25))
  cowincr25<-cbind(iterations,cowincr25)
  cowincr25['individual'] <- m
  cowincr25$algo <- 'ISAEM25'
  cowincr25 <- cowincr25[-1,]
  cowincr25$iterations = seq(1, end, by=1)
  final_25 <- rbind(final_25,cowincr25[0:end,])

  options_cowincr50<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
  cowincr50<-data.frame(saemix_incremental(saemix.model,saemix.data,options_cowincr50))
  cowincr50<-cbind(iterations,cowincr50)
  cowincr50['individual'] <- m
  cowincr50$algo <- 'ISAEM50'
  cow_scaled50 <- cowincr50[-1,]
  cow_scaled50$iterations = seq(1, 2*end, by=2)
  cow_scaled50 <- cow_scaled50[rep(seq_len(nrow(cow_scaled50)), each=2),]
  final_50 <- rbind(final_50,cow_scaled50[0:end,])
}


a <- graphConvMC_diffz(final_rwm[,c(1,2,9)],final_50[,c(1,2,9)],final_25[,c(1,2,9)])
b <- graphConvMC_diffw(final_rwm[,c(1,6,9)],final_50[,c(1,6,9)],final_25[,c(1,6,9)])

grid.arrange(a,b, ncol=2)
