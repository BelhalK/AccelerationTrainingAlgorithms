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
  # source('main_initialiseMainAlgo.R') 
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
source('initalgo.R') 


require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################


# Doc
oxboys.saemix<-read.table( "data/oxboys.saemix.tab",header=T,na=".")
oxboys.saemix_less <- oxboys.saemix[1:9,]
saemix.data<-saemixData(name.data=oxboys.saemix_less,header=TRUE,
  name.group=c("Subject"),name.predictors=c("age"),name.response=c("height"),
  units=list(x="yr",y="cm"))


growth.linear<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (2 columns, base and slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  base<-psi[id,1]
  slope<-psi[id,2]
  f<-base+slope*x
  return(f)
}
saemix.model<-saemixModel(model=growth.linear,description="Linear model",
  psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(1,0),covariance.model=matrix(c(1,1,1,1),ncol=2,byrow=TRUE), 
  error.model="constant")


indiv = 1
seed0 = 35644
replicate = 5
iter_mcmc = 10000
burn = 400


saemix.options_rwm<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,iter_mcmc,iter_mcmc,0))
saemix.options_linear<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))

#reference rwm
ref <- mcmc(saemix.model,saemix.data,saemix.options_rwm,iter_mcmc)
new<-mcmc(saemix.model,saemix.data,saemix.options_linear,iter_mcmc)




#expectations
expec_rwm <- ref$eta[[indiv]]
var_rwm <- ref$eta[[indiv]]
expec_rwm[,2:3] <- 0 
var_rwm[,2:3] <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_rwm<-list(seed=j+seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,iter_mcmc,iter_mcmc,0))
  post_rwm<-mcmc(saemix.model,saemix.data,saemix.options_rwm,iter_mcmc)$eta
  # print(post_rwm[[indiv]][44,2:3])
  post_rwm[[indiv]]['individual'] <- j
  expec_rwm[,2:3] <- expec_rwm[,2:3] + post_rwm[[indiv]][,2:3]
  var_rwm[,2] <- var_rwm[,2] + (post_rwm[[indiv]][,2])^2
  var_rwm[,3] <- var_rwm[,3] + (post_rwm[[indiv]][,3])^2
  
}
expec_rwm[,2:3] <- expec_rwm[,2:3]/replicate
var_rwm[,2:3] <- var_rwm[,2:3]/replicate

# graphConvMC_twokernels(expec_rwm,expec_rwm, title="Expectations")
# graphConvMC_twokernels(var_rwm,var_rwm, title="Variances")



expec_new <- new$eta[[indiv]]
var_new <- new$eta[[indiv]]
expec_new[,2:3] <- 0 
var_new[,2:3] <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_newkernel<-list(seed=j+seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))
  post_newkernel<-mcmc(saemix.model,saemix.data,saemix.options_newkernel,iter_mcmc)$eta
  post_newkernel[[indiv]]['individual'] <- j
  expec_new[,2:3] <- expec_new[,2:3] + post_newkernel[[indiv]][,2:3]
  var_new[,2] <- var_new[,2] + (post_newkernel[[indiv]][,2])^2
  var_new[,3] <- var_new[,3] + (post_newkernel[[indiv]][,3])^2
  
}
expec_new[,2:3] <- expec_new[,2:3]/replicate
var_new[,2:3] <- var_new[,2:3]/replicate


graphConvMC_twokernels(expec_rwm,expec_new, title="Expectations")
graphConvMC_twokernels(var_rwm,var_new, title="Variances")

#target is N(0,1)
#proposal is N(0,.01)
