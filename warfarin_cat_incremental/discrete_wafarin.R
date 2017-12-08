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
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat")
source('post_cat.R')
# source('main_cat.R')
source('main_cat_test.R')
source('main_estep_cat.R')
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


cat_data.saemix<-read.table("data/categorical1_data.txt",header=T,na=".")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("ID"),name.response=c("Y"),name.predictors=c("Y"), name.X=c("TIME"))


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


saemix.model<-saemixModel(model=cat_data.model,description="cat model",   
  psi0=matrix(c(0.5,0.4,0.3,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))),covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3, 
  byrow=TRUE),error.model="constant")



saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0))
saemix.foce<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc))


# post_rwm<-saemix_post_cat(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
# post_foce<-saemix_post_cat(saemix.model,saemix.data,saemix.foce)$post_newkernel


K1 = 100
K2 = 50

iterations = 1:(K1+K2+1)
gd_step = 0.01

#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(6,0,0,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE, map.range=c(0))
theo_ref<-data.frame(saemix_cat(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

graphConvMC_saem(theo_ref, title="new kernel")

#ref (map always)
options.cat<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(6,0,0,5),nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(3))
cat_saem<-data.frame(saemix_cat(saemix.model,saemix.data,options.cat))
cat_saem <- cbind(iterations, cat_saem)

graphConvMC_saem(cat_saem, title="new kernel")
graphConvMC2_saem(theo_ref,cat_saem, title="new kernel")


index = 1
graphConvMC_twokernels(post_rwm[[index]],post_rwm[[index]], title="rwm vs foce")
graphConvMC_twokernels(post_rwm[[index]],post_foce[[index]], title="rwm vs foce")


final_rwm <- post_rwm[[1]]
for (i in 2:length(post_rwm)) {
  final_rwm <- rbind(final_rwm, post_rwm[[i]])
}


final_foce <- post_foce[[1]]
for (i in 2:length(post_foce)) {
  final_foce <- rbind(final_foce, post_foce[[i]])
}



graphConvMC_twokernels(final_rwm,final_rwm, title="EM")
graphConvMC_twokernels(final_rwm,final_foce, title="EM")


#Autocorrelation
rwm.obj <- as.mcmc(post_rwm[[1]])
corr_rwm <- autocorr(rwm.obj[,2])
autocorr.plot(rwm.obj[,2])

foce.obj <- as.mcmc(post_foce[[1]])
corr_foce <- autocorr(foce.obj[,2])
autocorr.plot(foce.obj[,2])


#MSJD
mssd(post_rwm[[index]][,2])
mssd(post_foce[[index]][,2])



