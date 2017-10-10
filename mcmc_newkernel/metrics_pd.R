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
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mcmc_newkernel")
source('mcmc.R')
source('mcmc_mix.R')
source('mcmc_sum.R')



require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################



PD1.saemix<-read.table( "data/PD1.saemix.tab",header=T,na=".")
# PD1.saemix <- subset(PD1.saemix, dose!="90")

PD1.saemix_less <- PD1.saemix[295:297,]
PD2.saemix<-read.table( "data/PD2.saemix.tab",header=T,na=".")
saemix.data<-saemixData(name.data=PD1.saemix_less,header=TRUE,name.group=c("subject"),
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

indiv = 1
seed0 = 35644
replicate = 50
iter_mcmc = 20000
burn = 400


saemix.options_rwm<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,iter_mcmc,iter_mcmc,0))
saemix.options_linear<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))

#reference rwm
ref <- mcmc(saemix.model,saemix.data,saemix.options_rwm,iter_mcmc)
new<-mcmc(saemix.model,saemix.data,saemix.options_linear,iter_mcmc)

# #mix map and rwm
# saemix.options_mix<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,1,1,iter_mcmc))
# new_mix<-mcmc_mix(saemix.model,saemix.data,saemix.options_mix,iter_mcmc)

# #Sum of two proposals
# saemix.options_sum<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))
# new_sum<-mcmc_sum(saemix.model,saemix.data,saemix.options_sum,iter_mcmc)


graphConvMC_twokernels(ref$eta[[indiv]],ref$eta[[indiv]], title="eta")
graphConvMC_twokernels(new$eta[[indiv]],ref$eta[[indiv]], title="eta")
graphConvMC_twokernels(new$densy[[indiv]],ref$densy[[indiv]], title="Uy")
graphConvMC_twokernels(new$denseta[[indiv]],ref$denseta[[indiv]], title="Ueta")

graphConvMC_twokernels(new$denseta[[indiv]]+new$densy[[indiv]],ref$denseta[[indiv]]+ref$densy[[indiv]], title="Ueta")

graphConvMC_twokernels(new_sum$densy[[indiv]],ref$densy[[indiv]], title="Uy")




graphConvMC_twokernels(new$densy[[indiv]],ref$densy[[indiv]], title="Uy")
graphConvMC_twokernels(new_mix$densy[[indiv]],new$densy[[indiv]], title="Uy")
graphConvMC_twokernels(new_mix$densy[[indiv]],ref$densy[[indiv]], title="Uy")
graphConvMC_twokernels(new_mix$eta[[indiv]],ref$eta[[indiv]], title="eta")
graphConvMC_twokernels(new_mix$eta[[indiv]],new$eta[[indiv]], title="eta")

pack1 <- 100:230
pack2 <- 450:580
graphConvMC_twokernels(ref$densy[[indiv]][pack1,],ref$densy[[indiv]][pack2,], title="Uy")
graphConvMC_twokernels(ref$eta[[indiv]][pack1,],ref$eta[[indiv]][pack2,], title="eta")

graphConvMC_twokernels(new$densy[[indiv]][pack1,],new$densy[[indiv]][pack2,], title="Uy")
graphConvMC_twokernels(new$eta[[indiv]][pack1,],new$eta[[indiv]][pack2,], title="eta")

graphConvMC_twokernels(new_mix$densy[[indiv]][pack1,],new_mix$densy[[indiv]][pack2,], title="Uy")
graphConvMC_twokernels(new_mix$eta[[indiv]][pack1,],new_mix$eta[[indiv]][pack2,], title="eta")






#expectations
expec_rwm <- ref$eta[[indiv]]
var_rwm <- ref$eta[[indiv]]
expec_rwm[,2:4] <- 0 
var_rwm[,2:4] <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_rwm<-list(seed=j+seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0))
  post_rwm<-mcmc(saemix.model,saemix.data,saemix.options_rwm,iter_mcmc)$eta
  # print(post_rwm[[indiv]][44,2:4])
  post_rwm[[indiv]]['individual'] <- j
  expec_rwm[,2:4] <- expec_rwm[,2:4] + post_rwm[[indiv]][,2:4]
  var_rwm[,2] <- var_rwm[,2] + (post_rwm[[indiv]][,2])^2
  var_rwm[,3] <- var_rwm[,3] + (post_rwm[[indiv]][,3])^2
  var_rwm[,4] <- var_rwm[,4] + (post_rwm[[indiv]][,4])^2
}
expec_rwm[,2:4] <- expec_rwm[,2:4]/replicate
var_rwm[,2:4] <- var_rwm[,2:4]/replicate

# graphConvMC_twokernels(expec_rwm,expec_rwm, title="Expectations")
# graphConvMC_twokernels(var_rwm,var_rwm, title="Variances")



expec_new <- new$eta[[indiv]]
var_new <- new$eta[[indiv]]
expec_new[,2:4] <- 0 
var_new[,2:4] <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_newkernel<-list(seed=j+seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc))
  post_newkernel<-mcmc(saemix.model,saemix.data,saemix.options_newkernel,iter_mcmc)$eta
  post_newkernel[[indiv]]['individual'] <- j
  expec_new[,2:4] <- expec_new[,2:4] + post_newkernel[[indiv]][,2:4]
  var_new[,2] <- var_new[,2] + (post_newkernel[[indiv]][,2])^2
  var_new[,3] <- var_new[,3] + (post_newkernel[[indiv]][,3])^2
  var_new[,4] <- var_new[,4] + (post_newkernel[[indiv]][,4])^2
}
expec_new[,2:4] <- expec_new[,2:4]/replicate
var_new[,2:4] <- var_new[,2:4]/replicate


graphConvMC_twokernels(expec_new,expec_rwm, title="Expectations")
graphConvMC_twokernels(var_rwm,var_new, title="Variances")

#target is N(0,1)
#proposal is N(0,.01)
