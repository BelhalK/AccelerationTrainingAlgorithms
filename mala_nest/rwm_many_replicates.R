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
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mala_nest")
source("mixtureFunctions.R")
source('mala_main.R')
source('main_estep_mala.R')

library("mlxR")
require(ggplot2)
require(gridExtra)
require(reshape2)


#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
iter_mcmc = 800
replicate = 20
seed0 = 39546
indiv=4
burn = 500
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
  saemix.options_rwm<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0,0))
  post_rwm<-saemix_mala(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
  post_rwm[[indiv]]['individual'] <- j
  final_rwm <- rbind(final_rwm,post_rwm[[indiv]][-(1:burn),])
}


names(final_rwm)[1]<-paste("time")
names(final_rwm)[5]<-paste("id")
final_rwm <- final_rwm[c(5,1,2)]
prctilemlx(final_rwm[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1)
# graphConvMC_new(final_rwm, title="replicates")


final_mala <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_mala<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc,0))
  post_mala<-saemix_mala(saemix.model,saemix.data,saemix.options_mala)$post_mala
  post_mala[[indiv]]['individual'] <- j
  final_mala <- rbind(final_mala,post_mala[[indiv]][-(1:burn),])
}


names(final_mala)[1]<-paste("time")
names(final_mala)[5]<-paste("id")
final_mala <- final_mala[c(5,1,2)]
prctilemlx(final_mala[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1)
#ALl individual posteriors
# graphConvMC_new(final_mala, title="replicates")


final_nest <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_mala<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,iter_mcmc))
  post_nest<-saemix_mala(saemix.model,saemix.data,saemix.options_mala)$post_vb
  post_nest[[indiv]]['individual'] <- j
  final_nest <- rbind(final_nest,post_nest[[indiv]][-(1:burn),])
}


names(final_nest)[1]<-paste("time")
names(final_nest)[5]<-paste("id")
final_nest <- final_nest[c(5,1,2)]
prctilemlx(final_nest[-1,],band = list(number = 8, level = 80)) + ylim(-3,-1)
#ALl individual posteriors
# graphConvMC_new(final_nest, title="replicates")

prctilemlx(final_rwm[-1,],band = list(number = 8, level = 80))
prctilemlx(final_mala[-1,],band = list(number = 8, level = 80))
prctilemlx(final_nest[-1,],band = list(number = 8, level = 80))


final_rwm['algo'] <- 'rwm'
final_mala['algo'] <- 'mala'
final_nest['algo'] <- 'nest'
final <- rbind(final_rwm,final_mala)


labels <- c("rwm","mala")
prctilemlx(final) + theme(legend.position = "none")







