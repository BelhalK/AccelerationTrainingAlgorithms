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
source('mala_main.R')
source('main_estep_mala.R')
source("mixtureFunctions.R")


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

# Doc
theo.saemix<-read.table( "data/yield.saemix.tab",header=T,na=".")
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("site"),name.predictors=c("dose"),name.response=c("yield"),units=list(x="hr",y="mg/L",covariates=c("soil.nitrogen")), name.X="Time")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

yield.LP<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, ymax, xmax, slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  ymax<-psi[id,1]
  xmax<-psi[id,2]
  slope<-psi[id,3]
  f<-ymax+slope*(x-xmax)
#  cat(length(f),"  ",length(ymax),"\n")
  f[x>xmax]<-ymax[x>xmax]
  return(f)
}
saemix.model<-saemixModel(model=yield.LP,description="Linear plus plateau model",   
  psi0=matrix(c(8,100,0.2,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("Ymax","Xmax","slope"))),covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")



saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0,0,0,0))
saemix.options_mala<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc,0,0,0),sigma.val = 0.0001,gamma.val=0.001)
saemix.options_nest<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,0,iter_mcmc,0,0),sigma.val = 0.0001,gamma.val=0.001,memory=0.02)


post_rwm<-saemix_mala(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
post_mala<-saemix_mala(saemix.model,saemix.data,saemix.options_mala)$post_mala
post_nest<-saemix_mala(saemix.model,saemix.data,saemix.options_nest)$post_vb

index = 20
graphConvMC_twokernels(post_rwm[[index]],post_mala[[index]], title="RWM vs MALA")
graphConvMC_twokernels(post_mala[[index]],post_nest[[index]], title="MALA vs NEST")



names(post_rwm[[index]])[2]<-paste("ka")
names(post_rwm[[index]])[3]<-paste("V")
names(post_rwm[[index]])[4]<-paste("Cl")

names(post_mala[[index]])[2]<-paste("ka")
names(post_mala[[index]])[3]<-paste("V")
names(post_mala[[index]])[4]<-paste("Cl")

graphConvMC_twokernels(post_rwm[[index]][,-c(3,4)],post_mala[[index]][,-c(3,4)])


graphConvMC_twokernels(post_rwm[[index]],post_mala[[index]], title="RWM vs MALA")
graphConvMC_twokernels(post_mala[[index]],post_nest[[index]], title="MALA vs NEST")
graphConvMC_threekernels(post_rwm[[index]],post_mala[[index]],post_nest[[index]], title="EM")

final_rwm <- post_rwm[[1]]
for (i in 2:length(post_rwm)) {
  final_rwm <- rbind(final_rwm, post_rwm[[i]])
}


final_mala <- post_mala[[1]]
for (i in 2:length(post_mala)) {
  final_mala <- rbind(final_mala, post_mala[[i]])
}

final_nest <- post_nest[[1]]
for (i in 2:length(post_nest)) {
  final_nest <- rbind(final_nest, post_nest[[i]])
}


graphConvMC_threekernels(final_rwm,final_mala,final_nest, title="EM")


#ALl individual posteriors
graphConvMC_new(final_rwm, title="RWM")
graphConvMC_new(final_mala, title="MALA")
graphConvMC_new(final_nest, title="Nest")

#first individual posteriors
graphConvMC_new(post_rwm[[1]], title="EM")

graphConvMC_twokernels(final_rwm,final_mala, title="EM")
graphConvMC_threekernels(final_rwm,final_mala,final_nest, title="EM")
graphConvMC_twokernels(post_rwm[[1]],post_mala[[1]], title="EM")





# theo.onlypop<-saemix(saemix.model,saemix.data,saemix.options)

# saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
# plot(saemix.fit,plot.type="individual")

