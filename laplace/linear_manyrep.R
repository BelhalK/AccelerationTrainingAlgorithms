#library(rstan)
setwd("/Users/karimimohammedbelhal/Desktop/variationalBayes/mcmc_R_isolate/Dir2")
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_cov.R') 
  source('func_distcond.R') 
  source('func_FIM.R') 
  source('func_ggplot2.R') 
  source('func_plots.R') 
  source('func_simulations.R') 
  source('ggplot2_global.R') 
  # source('KL.R') 
  #source('vi.R') 
  source('global.R')
  source('main.R')
  source('mcmc_main.R') 
  source('main_estep.R')
  source('main_estep_mcmc.R') 
  source('main_estep_morekernels.R') 
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('plots_ggplot2.R') 
  source('saemix-package.R') 
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/laplace")
source('laplace_main.R')
source('main_estep_laplace.R')
source("mixtureFunctions.R")

library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
iter_mcmc = 200
replicate = 5
seed0 = 39546
indiv=2
burn = 50




theo.saemix<-read.table("data/linear_matlab2.txt",header=TRUE,na=".",sep=",")
# theo.saemix<-read.table("data/linear_matlab.txt",header=TRUE,na=".",sep=",")
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Time"),name.response=c("y"),name.X="Time")

model1cpt<-function(psi,id,xidep) { 
  tim<-xidep[,1]  
  d<-psi[id,1]
  b<-psi[id,2]

  ypred<-d*tim+b
  return(ypred)
}
# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(5,5),ncol=2,byrow=TRUE, dimnames=list(NULL, c("d","b"))),transform.par=c(0,0))

saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0,0,0,0))
saemix.laplace<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc,0,0,0))
# saemix.fo<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,iter_mcmc,0,0))
saemix.fo2<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,0,iter_mcmc,0))
saemix.foce<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,0,0,iter_mcmc))


post_rwm<-saemix_laplace(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
post_foce<-saemix_laplace(saemix.model,saemix.data,saemix.foce)$post_newkernel
post_laplace<-saemix_laplace(saemix.model,saemix.data,saemix.laplace)$post_newkernel
# post_fo<-saemix_laplace(saemix.model,saemix.data,saemix.fo)$post_newkernel
post_fo2<-saemix_laplace(saemix.model,saemix.data,saemix.fo2)$post_newkernel


final_rwm <- 0
for (j in 1:replicate){
  print(j)
  saemix.rwm<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0,0,0,0))
  post_rwm<-saemix_laplace(saemix.model,saemix.data,saemix.rwm)$post_rwm
  post_rwm[[indiv]]['individual'] <- j
  final_rwm <- rbind(final_rwm,post_rwm[[indiv]][-1,])
}


names(final_rwm)[1]<-paste("time")
names(final_rwm)[4]<-paste("id")
final_rwm <- final_rwm[c(4,1,2)]
# prctilemlx(final_rwm[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1) + ggtitle("RWM")

#burn
rwm_burn <- final_rwm[final_rwm[,2]>burn,]


final_foce <- 0
for (j in 1:replicate){
  print(j)
  saemix.foce<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,0,0,iter_mcmc))
  post_foce<-saemix_laplace(saemix.model,saemix.data,saemix.foce)$post_newkernel
  post_foce[[indiv]]['individual'] <- j
  final_foce <- rbind(final_foce,post_foce[[indiv]][-1,])
}


names(final_foce)[1]<-paste("time")
names(final_foce)[4]<-paste("id")
final_foce <- final_foce[c(4,1,2)]
# prctilemlx(final_foce[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1) + ggtitle("foce")

#burn
foce_burn <- final_foce[final_foce[,2]>burn,]



final_fo2 <- 0
for (j in 1:replicate){
  print(j)
  saemix.fo2<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,0,iter_mcmc,0))
  post_fo2<-saemix_laplace(saemix.model,saemix.data,saemix.fo2)$post_newkernel
  post_fo2[[indiv]]['individual'] <- j
  final_fo2 <- rbind(final_fo2,post_fo2[[indiv]][-1,])
}


names(final_fo2)[1]<-paste("time")
names(final_fo2)[4]<-paste("id")
final_fo2 <- final_fo2[c(4,1,2)]
# prctilemlx(final_fo2[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1) + ggtitle("fo2")

#burn
fo2_burn <- final_fo2[final_fo2[,2]>burn,]



final_laplace <- 0
for (j in 1:replicate){
  print(j)
  saemix.laplace<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc,0,0,0))
  post_laplace<-saemix_laplace(saemix.model,saemix.data,saemix.laplace)$post_newkernel
  post_laplace[[indiv]]['individual'] <- j
  final_laplace <- rbind(final_laplace,post_laplace[[indiv]][-1,])
}


names(final_laplace)[1]<-paste("time")
names(final_laplace)[4]<-paste("id")
final_laplace <- final_laplace[c(4,1,2)]
# prctilemlx(final_laplace[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1) + ggtitle("laplace")

#burn
laplace_burn <- final_laplace[final_laplace[,2]>burn,]



rwm_burn['group'] <- 1
foce_burn['group'] <- 2
foce_burn$id <- foce_burn$id +1
fo2_burn['group'] <- 3
fo2_burn$id <- fo2_burn$id +2
laplace_burn['group'] <- 4
laplace_burn$id <- laplace_burn$id +3

final <- 0
final <- rbind(rwm_burn,foce_burn, fo2_burn,laplace_burn)


labels <- c("rwm","foce","fo2","laplace")
final <- final[c(1,4,2,3)]
prctilemlx(final, band = list(number = 2, level = 80),group='group', label = labels) + theme(legend.position = "none")




#Autocorrelation
rwm.obj <- as.mcmc(post_rwm[[1]])
corr_rwm <- autocorr(rwm.obj[,2])
autocorr.plot(rwm.obj[,2],main="rwm")

foce.obj <- as.mcmc(post_foce[[1]])
corr_foce <- autocorr(foce.obj[,2])
autocorr.plot(foce.obj[,2])

fo2.obj <- as.mcmc(post_fo2[[1]])
corr_fo2 <- autocorr(fo2.obj[,2])
autocorr.plot(fo2.obj[,2])

laplace.obj <- as.mcmc(post_laplace[[1]])
corr_laplace <- autocorr(laplace.obj[,2])
autocorr.plot(laplace.obj[,2])



#MSJD
print(paste0("msjd rwm: ", mssd(rwm_burn[,3])))

mssd(foce_burn[,3])
mssd(fo2_burn[,3])
mssd(laplace_burn[,3])



