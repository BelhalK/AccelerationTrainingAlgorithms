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

# Doc
data(yield.saemix)
yield.saemix_less <- yield.saemix[1:6,]
saemix.data<-saemixData(name.data=yield.saemix_less,header=TRUE,name.group=c("site"),
  name.predictors=c("dose"),name.response=c("yield"),
  name.covariates=c("soil.nitrogen"),units=list(x="kg/ha",y="t/ha",
  covariates=c("kg/ha")))

yield.LP<-function(psi,id,xidep) {
  x<-xidep[,1]
  ymax<-psi[id,1]
  xmax<-psi[id,2]
  slope<-psi[id,3]
  f<-ymax+slope*(x-xmax)
  f[x>xmax]<-ymax[x>xmax]
  return(f)
}

saemix.model<-saemixModel(model=yield.LP,description="Linear plus plateau model",   
  psi0=matrix(c(8,100,0.2,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("Ymax","Xmax","slope"))),covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")


indiv = 1
seed0 = 35644
replicate = 5
iter_mcmc = 1000
burn = 400


saemix.options_rwm<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,iter_mcmc,iter_mcmc,0))
saemix.options_linear<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))

ref <- mcmc(saemix.model,saemix.data,saemix.options_rwm,iter_mcmc)
new<-mcmc(saemix.model,saemix.data,saemix.options_linear,iter_mcmc)

saemix.options_mix<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,1,1,iter_mcmc))
new_mix<-mcmc_mix(saemix.model,saemix.data,saemix.options_mix,iter_mcmc)


saemix.options_sum<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))
new_sum<-mcmc_sum(saemix.model,saemix.data,saemix.options_sum,iter_mcmc)


graphConvMC_twokernels(new$eta[[indiv]],ref$eta[[indiv]], title="eta")
graphConvMC_twokernels(new$densy[[indiv]],ref$densy[[indiv]], title="Uy")
graphConvMC_twokernels(new$denseta[[indiv]],ref$denseta[[indiv]], title="Ueta")

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

