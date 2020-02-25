setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/PackageContributions/FINALPACKAGES/saemixextension/R")
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
setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/PackageContributions/FINALPACKAGES/saemixextension/demos")





###WARFA
smx.ord <- read.table("data/categorical1_data.txt",header=T)
data.ord<-saemixData(name.data=smx.ord,header=TRUE,sep=" ",na=NA, 
name.group=c("ID"),name.response=c("Y"),name.predictors=c("TIME","Y"), 
name.X=c("TIME"))

ordinal.model<-function(psi,id,xidep) {
   T<-xidep[,1]
   y<-xidep[,2]
   alp1<-psi[id,1]
   alp2<-psi[id,2]
   alp3<-psi[id,3]
   logit1<-alp1
   logit2<-alp1+alp2
   logit3<-alp1+alp2+alp3
   pge1<-exp(logit1)/(1+exp(logit1))
   pge2<-exp(logit2)/(1+exp(logit2))
   pge3<-exp(logit3)/(1+exp(logit3))
   logpdf<-rep(0,length(T))
    P.obs = (y==0)*pge1+(y==1)*(pge2 - pge1)+(y==2)*(pge3 - pge2)+(y==3)*(1 - pge3)
    logpdf <- log(P.obs)

   # logpdf[y==0]<-log(pge1)
   # logpdf[y==1]<-log(pge2 - pge1)
   # logpdf[y==2]<-log(pge3 - pge2)
   # logpdf[y==3]<-log(1 - pge3)

   return(logpdf)
}


model.ord<-saemixModel(model=ordinal.model,description="Ordinal categorical 
model",type="likelihood",
 
psi0=matrix(c(3,1,0.5),ncol=3,byrow=TRUE,dimnames=list(NULL,c("alp1","alp2","alp3"))),
 
transform.par=c(0,1,1),covariance.model=matrix(c(1,rep(0,8)),ncol=3))

options.ord<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), 
nbiter.saemix = c(200,100), nbiter.sa=0, 
displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)

##RUNS
K1 = 200
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)
catordinal<-saemix(model.ord,data.ord,options.ord)


test <- predict(catordinal)
test <- saemix.predict(catordinal)
