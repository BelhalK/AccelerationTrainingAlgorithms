#library(rstan)
setwd("/home/belhal.karimi/Desktop/Belhal/Dir2")
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
setwd("/home/belhal.karimi/Desktop/Belhal/mcmc_newkernel")
source('mcmc.R')



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
iter_mcmc = 1000000
burn = 400


saemix.options_rwm<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,iter_mcmc,iter_mcmc,0))
saemix.options_linear<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))

# ref <- mcmc(saemix.model,saemix.data,saemix.options_rwm,iter_mcmc)
new<-mcmc(saemix.model,saemix.data,saemix.options_linear,iter_mcmc)

  
# eta <- graphConvMC_twokernels(new$eta[[indiv]],ref$eta[[indiv]], title="eta")
# ggsave(plot = eta, file = paste("eta.pdf"))

pack1 <- 10000:20000
pack2 <- 100000:110000
pack3 <- 200000:210000
pack4 <- 300000:310000
pack5 <- 400000:410000
pack6 <- 500000:510000
pack7 <- 600000:610000
pack8 <- 700000:710000
pack9 <- 800000:810000
pack10 <- 900000:910000


Uy1 <- graphConvMC_twokernels(new$densy[[indiv]][pack1,],new$densy[[indiv]][pack1,], title="Uy")
ggsave(plot = Uy1, file = paste("densy1.pdf"))
Uy2 <- graphConvMC_twokernels(new$densy[[indiv]][pack2,],new$densy[[indiv]][pack2,], title="Uy")
ggsave(plot = Uy2, file = paste("densy2.pdf"))
Uy3 <- graphConvMC_twokernels(new$densy[[indiv]][pack3,],new$densy[[indiv]][pack3,], title="Uy")
ggsave(plot = Uy3, file = paste("densy3.pdf"))
Uy4 <- graphConvMC_twokernels(new$densy[[indiv]][pack4,],new$densy[[indiv]][pack4,], title="Uy")
ggsave(plot = Uy4, file = paste("densy4.pdf"))
Uy5 <- graphConvMC_twokernels(new$densy[[indiv]][pack5,],new$densy[[indiv]][pack5,], title="Uy")
ggsave(plot = Uy5, file = paste("densy5.pdf"))
Uy6 <- graphConvMC_twokernels(new$densy[[indiv]][pack6,],new$densy[[indiv]][pack6,], title="Uy")
ggsave(plot = Uy6, file = paste("densy6.pdf"))
Uy7 <- graphConvMC_twokernels(new$densy[[indiv]][pack7,],new$densy[[indiv]][pack7,], title="Uy")
ggsave(plot = Uy7, file = paste("densy7.pdf"))
Uy8 <- graphConvMC_twokernels(new$densy[[indiv]][pack8,],new$densy[[indiv]][pack8,], title="Uy")
ggsave(plot = Uy8, file = paste("densy8.pdf"))
Uy9 <- graphConvMC_twokernels(new$densy[[indiv]][pack9,],new$densy[[indiv]][pack9,], title="Uy")
ggsave(plot = Uy9, file = paste("densy9.pdf"))
Uy10 <- graphConvMC_twokernels(new$densy[[indiv]][pack10,],new$densy[[indiv]][pack10,], title="Uy")
ggsave(plot = Uy10, file = paste("densy10.pdf"))


# Ueta <- graphConvMC_twokernels(new$denseta[[indiv]],ref$denseta[[indiv]], title="Ueta")
# ggsave(plot = Ueta, file = paste("denseta.pdf"))



