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
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem")
source('newkernel_main.R')
source('main_new.R')
source('main_estep_new.R')
source('main_gd.R')
source('main_estep_gd.R')
source('main_gd_mix.R')
source('main_estep_gd_mix.R')
source('main_estep_mix.R')
source('main_estep_newkernel.R')
source("mixtureFunctions.R")

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/laplace_incremental_saem")
source('main_laplace_incremental.R')
source('main_estep_laplace_incremental.R')




PD1.saemix<-read.table( "PD1.saemix.tab",header=T,na=".")
PD1.saemix <- subset(PD1.saemix, dose!="90")


PD2.saemix<-read.table( "PD2.saemix.tab",header=T,na=".")
saemix.data1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))
saemix.data2<-saemixData(name.data=PD2.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))


modelemax<-function(psi,id,xidep) {
# input:
# psi : matrix of parameters (3 columns, E0, Emax, EC50)
# id : vector of indices
# xidep : dependent variables (same nb of rows as length of id)
# returns:
# a vector of predictions of length equal to length of id
dose<-xidep[,1]
e0<-psi[id,1]
emax<-psi[id,2]
e50<-psi[id,3]
f<-e0+emax*dose/(e50+dose)
return(f)
}


saemix.model<-saemixModel(model=modelemax,description="Emax model",
psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")

K1 = 150
K2 = 30
iteration = 1:(K1+K2+1)




#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0)
theo_ref<-data.frame(saemix(saemix.model,saemix.data,options))
theo_ref <- cbind(iteration, theo_ref)

#ref (map always)
# options.new<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5),nbiter.saemix = c(K1,K2))
# theo_new_ref<-data.frame(saemix_new(saemix.model,saemix.data,options.new))
# theo_new_ref <- cbind(iteration, theo_new_ref)


#FOCE
options.mix<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,4,0),nbiter.saemix = c(K1,K2),step.gd=gd_step,map.range=c(3))
theo_mix<-data.frame(saemix_gd_mix(saemix.model,saemix.data,options.mix))
theo_mix <- cbind(iteration, theo_mix)

#FOCE
options.lap<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,4), nbiter.saemix = c(K1,K2),step.gd=gd_step,map.range=c(3),nb.replacement=100)
theo_lap <- data.frame(saemix_laplace_incremental(saemix.model,saemix.data,options.lap))
theo_lap <- cbind(iteration, theo_lap)




theo_ref$algo <- 'rwm'
theo_incremental$algo <- 'ISAEM'

theo_ref_scaled <- theo_ref[rep(seq_len(nrow(theo_ref)), each=2),]
theo_ref_scaled$iteration = 1:(2*(K1+K2+1))


comparison <- 0
# comparison <- rbind(theo_ref,theo_incremental)
comparison <- rbind(theo_ref_scaled[iteration,],theo_incremental)

var <- melt(comparison, id.var = c('iteration','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)


