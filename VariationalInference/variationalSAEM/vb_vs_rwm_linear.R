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
setwd("/Users/karimimohammedbelhal/Desktop/git_phd/saem/VariationalInference/variationalSAEM")
source('vb_main.R')
source('main_estep_vb.R')

require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
iter_mcmc = 10

# Doc
theo.saemix<-read.table("data/linear_matlab.txt",header=TRUE,na=".",sep=",")
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Time"),name.response=c("y"),name.X="Time")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

model1cpt<-function(psi,id,xidep) { 
	tim<-xidep[,1]  
	d<-psi[id,1]
	b<-psi[id,2]

	ypred<-d*tim+b
	return(ypred)
}
# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(5,5),ncol=2,byrow=TRUE, dimnames=list(NULL, c("d","b"))),transform.par=c(0,0))

saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0,0))
saemix.options_linear<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc,0))
saemix.options_vb<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,0,iter_mcmc))

post_rwm<-saemix_vb(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
post_vb_linear<-saemix_vb(saemix.model,saemix.data,saemix.options_linear)$post_vb_linear
post_vb<-saemix_vb(saemix.model,saemix.data,saemix.options_vb)$post_vb


final_rwm <- post_rwm[[1]]
for (i in 2:length(post_rwm)) {
  final_rwm <- rbind(final_rwm, post_rwm[[i]])
}

final_vb <- post_vb[[1]]
for (i in 2:length(post_vb)) {
  final_vb <- rbind(final_vb, post_vb[[i]])
}

final_vb_linear <- post_vb_linear[[1]]
for (i in 2:length(post_vb_linear)) {
  final_vb_linear <- rbind(final_vb, post_vb_linear[[i]])
}

#ALl individual posteriors
graphConvMC_new(final_rwm, title="RWM")
graphConvMC_new(final_vb, title="Variational Bayes")
graphConvMC_new(final_vb_linear, title="Variational Bayes")
graphConvMC_twokernels(final_rwm,final_vb, title="EM")
#first individual posteriors
graphConvMC_new(post_rwm[[1]], title="EM")
graphConvMC_twokernels(post_rwm[[1]],post_vb[[1]], title="EM")





# theo.onlypop<-saemix(saemix.model,saemix.data,saemix.options)

# saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
# plot(saemix.fit,plot.type="individual")

