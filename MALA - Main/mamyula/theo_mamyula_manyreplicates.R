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
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mamyula")
source("mixtureFunctions.R")
source('mala_main.R')
source('main_estep_mala.R')

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
iter_mcmc = 1000
replicate = 2
seed0 = 39546
indiv=4
burn = 600
# Doc
theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
l <- c(4.02,4.4,4.53,4.4,5.86,4,4.95,4.53,3.1,5.5,4.92,5.3)
for (i in 1:12){
  theo.saemix[(i*10-9):(i*10),'Dose'] = l[i]
}
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

model1cpt<-function(psi,id,xidep) { 
	dose<-xidep[,1]
	tim<-xidep[,2]  
	ka<-psi[id,1]
	V<-psi[id,2]
	CL<-psi[id,3]
	k<-CL/V
	ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
	return(ypred)
}

# Default model, no covariate
# saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))
saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(10,10,1.05),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))

final_rwm <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_rwm<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0,0,0,0))
  post_rwm<-saemix_mala(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
  post_rwm[[indiv]]['individual'] <- j
  final_rwm <- rbind(final_rwm,post_rwm[[indiv]][-1,])
}


names(final_rwm)[1]<-paste("time")
names(final_rwm)[5]<-paste("id")
final_rwm <- final_rwm[c(5,1,2)]
# prctilemlx(final_rwm[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1) + ggtitle("RWM")

#burn
rwm_burn <- final_rwm[final_rwm[,2]>burn,]
prctilemlx(rwm_burn[-1,-4],band = list(number = 2, level = 80)) + ylim(-3,-1) + ggtitle("RWM")



final_mala <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_mala<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc,0,0,0),sigma.val=0.01,gamma.val=0.01)
  post_mala<-saemix_mala(saemix.model,saemix.data,saemix.options_mala)$post_mala
  post_mala[[indiv]]['individual'] <- j
  final_mala <- rbind(final_mala,post_mala[[indiv]][-1,])
}


names(final_mala)[1]<-paste("time")
names(final_mala)[5]<-paste("id")
final_mala <- final_mala[c(5,1,2)]
# prctilemlx(final_mala[-1,],band = list(number = 4, level = 80)) + ylim(-3,-1) + ggtitle("MALA")

#burn
mala_burn <- final_mala[final_mala[,2]>burn,]
prctilemlx(mala_burn[-1,-4],band = list(number = 2, level = 80)) + ylim(-3,-1) + ggtitle("MALA")


final_mamyula <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_mamyula<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,iter_mcmc,0,0),sigma.val=0.1,gamma.val=0.01)
  post_mamyula<-saemix_mala(saemix.model,saemix.data,saemix.options_mamyula)$post_mala
  post_mamyula[[indiv]]['individual'] <- j
  final_mamyula <- rbind(final_mamyula,post_mamyula[[indiv]][-1,])
}



names(final_mamyula)[1]<-paste("time")
names(final_mamyula)[5]<-paste("id")
final_mamyula <- final_mamyula[c(5,1,2)]
# prctilemlx(final_mamyula[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1) + ggtitle("mamyulaerov")



#burn
mamyula_burn <- final_mamyula[final_mamyula[,2]>burn,]
prctilemlx(mamyula_burn[-1,],band = list(number = 2, level = 80)) + ylim(-3,-1) + ggtitle("mamyula")



#ALl individual posteriors
# graphConvMC_new(final_mamyula, title="replicates")


rwm_burn['group'] <- 1
mala_burn['group'] <- 2
mala_burn$id <- mala_burn$id +1
mamyula_burn['group'] <- 3
mamyula_burn$id <- mamyula_burn$id +2

final <- 0
final <- rbind(rwm_burn,mala_burn)
final <- rbind(rwm_burn,mala_burn, mamyula_burn)


labels <- c("rwm","mala","mamyula")
labels <- c("rwm","mala")
final <- final[c(1,4,2,3)]
prctilemlx(final, band = list(number = 2, level = 80),group='group', label = labels) + theme(legend.position = "none")


final_rwm <- final_rwm[,-c(4,5)]
final_mala <- final_mala[,-c(4,5)]
final_mamyula <- final_mamyula[,-c(4,5)]



#Autocorrelation
rwm.obj <- as.mcmc(post_rwm[[1]])
corr_rwm <- autocorr(rwm.obj[,2])
autocorr.plot(rwm.obj[,2])

mala.obj <- as.mcmc(post_mala[[1]])
corr_mala <- autocorr(mala.obj[,2])
autocorr.plot(mala.obj[,2])

mamyula.obj <- as.mcmc(post_mamyula[[1]])
corr_mamyula <- autocorr(mamyula.obj[,2])
autocorr.plot(mamyula.obj[,2])

#MSJD
mssd(rwm_burn[,3])
mssd(mala_burn[,3])
mssd(mamyula_burn[,3])
