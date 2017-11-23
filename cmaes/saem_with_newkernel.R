
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/cmaes/Dir")
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
  source('main_new.R')
  source('main_estep_new2.R')
  source('main_new_mix.R')
  source('main_estep_mix.R')
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/cmaes/")
source("mixtureFunctions.R")


library("rJava")
library("rCMA")
library("mlxR")
library(sgd)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)


#####################################################################################
# Theophylline

# Data - changing gender to M/F
theo.saemix<-read.table("theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")



# theo.saemix_less <- theo.saemix[1:120,]
theo.saemix_less <- theo.saemix
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
  ,covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3,byrow=TRUE),psi0=matrix(c(1.,20,0.5),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))

K1 = 20
K2 = 10
iterations = 1:(K1+K2+1)
gd_step = 0.01

options<-list(seed=395246,map=F,fim=F,ll.is=T,displayProgress=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,0),nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0,cma=0)
war_ref<-data.frame(saemix_new_mix(saemix.model,saemix.data,options))
war_ref <- cbind(iterations, war_ref)

options.cma<-list(seed=395246,map=F,fim=F,ll.is=T,displayProgress=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,0),nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0,cma=1)
war_cma<-data.frame(saemix_new_mix(saemix.model,saemix.data,options.cma))
war_cma <- cbind(iterations, war_cma)

war_ref[(K1+K2),]
war_cma[(K1+K2),]


graphConvMC_twokernels(war_ref,war_cma, title="new kernel")



options.new<-list(seed=395246,map=F,fim=F,ll.is=T,displayProgress=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),map.range=c(1:3),nbiter.burn =0)
theo_new_ref<-data.frame(saemix_new_mix(saemix.model,saemix.data,options.new))
theo_new_ref <- cbind(iterations, theo_new_ref)
theo_new_ref[end,]
graphConvMC_twokernels(theo_new_ref,theo_new_ref, title="new kernel")