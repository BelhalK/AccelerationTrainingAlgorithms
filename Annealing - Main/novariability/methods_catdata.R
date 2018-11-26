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

# library('rCMA')
# source("/Users/karimimohammedbelhal/Desktop/papers/iem_code/imcem_saemix/plots_se.R")


# library("rJava")
# library("rCMA")
# library("mlxR")
# library("psych")
# library("coda")
# library("Matrix")
# library(abind)
# require(ggplot2)
# require(gridExtra)
# require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
iter_mcmc = 200


# cat_data.saemix<-read.table("data/categorical1_data.txt",header=T,na=".")
# cat_data.saemix <- cat_data.saemix[1:692,]
# saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("ID"),name.response=c("Y"),name.predictors=c("Y"), name.X=c("TIME"))



# cat_data.saemix<-read.table("data/categorical1_data.txt",header=T,na=".")
# cat_data.saemix<-read.table("data/categorical1_data_less2.txt",header=T,na=".")
# cat_data.saemix<-read.table("data/categorical1_data_less2.txt",header=T,na=".")
# cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat1.csv", header=T, sep=",")
cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat2.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("y"), name.X=c("time"))
# saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("ID"),name.response=c("Y"),name.predictors=c("Y"), name.X=c("TIME"))



cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]

th1 <- psi[id,1]
th2 <- psi[id,2]
th3 <- psi[id,3]

P0 <- 1/(1+exp(-th1))
Pcum1 <- 1/(1+exp(-th1-th2))
Pcum2 <- 1/(1+exp(-th1-th2-th3))

P1 <- Pcum1 - P0
P2 <- Pcum2 - Pcum1
P3 <- 1 - Pcum2

P.obs = (level==0)*P0+(level==1)*P1+(level==2)*P2+(level==3)*P3

return(P.obs)
}


saemix.model<-saemixModel(model=cat_data.model,description="cat model",type="likelihood",   
  psi0=matrix(c(2,1,2),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))), 
  transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),omega.init=matrix(c(2,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")


saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0))
saemix.foce<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc))


# post_rwm<-saemix_post_cat(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
# post_foce<-saemix_post_cat(saemix.model,saemix.data,saemix.foce)$post_newkernel


K1 = 500
K2 = 100

iterations= 1:(K1+K2+1)
gd_step = 0.01
end = K1+K2
seed0 = 444

#RWM
replicate <- 3
final_optim <- 0
final_av <- 0
final_avnew <- 0
final_bayes <- 0
for (m in 1:replicate){
  print(m)
  print(m)
  # l = list(c(1,5,1,0,0,0),c(0.8,12,0.8,0,0,0),c(1.2,3,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  l = list(c(2,1,2),c(2,1,5),c(2,1,3))
  saemix.model<-saemixModel(model=cat_data.model,description="cat model", type="likelihood",   
  psi0=matrix(l[[m]],ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))), 
  transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE),omega.init=matrix(c(3/m,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")


  #No var
  ##### Optim (fmin search)
  options_cat_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
  cat_optim<-data.frame(saemix(saemix.model,saemix.data,options_cat_with))
  cat_optim <- cbind(iterations, cat_optim)
  cat_optim['individual'] <- m
  final_optim <- rbind(final_optim,cat_optim)
  ##### AV
  options_cat_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,0),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, map.range=c(0),av=1)
  cat_withav<-data.frame(saemix(saemix.model,saemix.data,options_cat_with))
  cat_withav <- cbind(iterations, cat_withav)
  cat_withav['individual'] <- m
  final_av <- rbind(final_av,cat_withav)

  ##### AV and new kernel
  options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:3), av=1)
  cat_newkernelav<-data.frame(saemix(saemix.model,saemix.data,options_newkernel))
  cat_newkernelav <- cbind(iterations, cat_newkernelav)
  cat_newkernelav['individual'] <- m
  final_avnew <- rbind(final_avnew,cat_newkernelav)

  ##### pseudo bayesian
  options_cat_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0,map.range=c(0))
  cat_bayes<-data.frame(saemix(saemix.model,saemix.data,options_cat_with))
  cat_bayes <- cbind(iterations, cat_bayes)
  cat_bayes['individual'] <- m
  final_bayes <- rbind(final_bayes,cat_bayes)
}

graphConvMC_diff4(final_optim,final_av,final_avnew,final_bayes, title="")

#black: optim
#blue: av
#red: av newkernel
#green: bayes

graphConvMC_diff(final_optim,final_bayes, title="")