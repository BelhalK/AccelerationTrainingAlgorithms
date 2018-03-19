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
  source("mixtureFunctions.R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel")
source('newkernel_main.R')
source('main_estep_newkernel.R')

library("mlxR")
library("psych")
library("coda")
library("Matrix")

require(ggplot2)
require(gridExtra)
require(reshape2)


#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
iter_mcmc = 700
replicate = 2
seed0 = 39546
indiv=4
burn = 300
# Doc
# Doc
data(theo.saemix)
theo.saemix_less <- theo.saemix[1:120,]
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
saemix.data<-saemixData(name.data=theo.saemix_less,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

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
# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption"
  ,psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))



final_rwm <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_rwm<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0))
  post_rwm<-saemix_newkernel(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
  post_rwm[[indiv]]['individual'] <- j
  final_rwm <- rbind(final_rwm,post_rwm[[indiv]][-1,])
}


names(final_rwm)[1]<-paste("time")
names(final_rwm)[5]<-paste("id")
final_rwm <- final_rwm[c(5,1,2)]
# prctilemlx(final_rwm[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1) + ggtitle("RWM")

#burn
rwm_burn <- final_rwm[final_rwm[,2]>burn,]
prctilemlx(rwm_burn[-1,],band = list(number = 4, level = 80)) + ylim(-3,-1) + ggtitle("RWM")



final_new <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_new<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc))
  post_new<-saemix_newkernel(saemix.model,saemix.data,saemix.options_new)$post_newkernel
  post_new[[indiv]]['individual'] <- j
  final_new <- rbind(final_new,post_new[[indiv]][-1,])
}


names(final_new)[1]<-paste("time")
names(final_new)[5]<-paste("id")
final_new <- final_new[c(5,1,2)]
# prctilemlx(final_new[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1) + ggtitle("RWM")

#burn
new_burn <- final_new[final_new[,2]>burn,]
prctilemlx(new_burn[-1,],band = list(number = 4, level = 80)) + ylim(-3,-1) + ggtitle("new")


# graphConvMC_twokernels(post_rwm[[1]],post_new[[1]], title="EM")

final_rwm <- post_rwm[[1]]
for (i in 2:length(post_rwm)) {
  final_rwm <- rbind(final_rwm, post_rwm[[i]])
}


final_new <- post_new[[1]]
for (i in 2:length(post_new)) {
  final_new <- rbind(final_new, post_new[[i]])
}


rwm.obj <- as.mcmc(post_rwm[[1]])
corr_rwm <- autocorr(rwm.obj[,2])
autocorr.plot(rwm.obj[,2]) + title("RWM SAEM Autocorrelation")

new.obj <- as.mcmc(post_new[[1]])
corr_new <- autocorr(new.obj[,2])
autocorr.plot(new.obj[,2]) + title("Laplace SAEM Autocorrelation")


#MSJD
mssd(rwm_burn[,3])
mssd(new_burn[,3])
