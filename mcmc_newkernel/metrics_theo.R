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
source('initalgo.R') 
source('mcmc_sum.R')
library("mlxR")

require(ggplot2)
require(gridExtra)
require(reshape2)



#####################################################################################

# Doc
data(theo.saemix)
theo.saemix_less <- theo.saemix[1:120,]
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
  ,psi0=matrix(c(1,80,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))


indiv = 1
seed0 = 35644
replicate = 20
iter_mcmc = 1000
burn = 400


saemix.options_rwm<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,iter_mcmc,iter_mcmc,0))
saemix.options_linear<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))

#reference rwm
ref <- mcmc(saemix.model,saemix.data,saemix.options_rwm,iter_mcmc)
new<-mcmc(saemix.model,saemix.data,saemix.options_linear,iter_mcmc)

saemix.options_linear<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))
new<-mcmc(saemix.model,saemix.data,saemix.options_linear,iter_mcmc)

saemix.options_linear<-list(seed=seed0+44,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))
new2<-mcmc(saemix.model,saemix.data,saemix.options_linear,iter_mcmc)
graphConvMC_twokernels(new$eta[[indiv]],new2$eta[[indiv]], title="eta")


graphConvMC_twokernels(ref$eta[[indiv]],ref$eta[[indiv]], title="eta")
graphConvMC_twokernels(new$eta[[indiv]],new$eta[[indiv]], title="eta")

graphConvMC_twokernels(ref$eta[[indiv]],new$eta[[indiv]], title="eta")

final_rwm <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_rwm<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,iter_mcmc,iter_mcmc,0))
  post_rwm<-mcmc(saemix.model,saemix.data,saemix.options_rwm)$eta[[indiv]]
  post_rwm['individual'] <- j
  final_rwm <- rbind(final_rwm,post_rwm)
}




names(final_rwm)[1]<-paste("time")
names(final_rwm)[5]<-paste("id")
final_rwm1 <- final_rwm[c(5,1,2)]
final_rwm2 <- final_rwm[c(5,1,3)]
final_rwm3 <- final_rwm[c(5,1,4)]


prctilemlx(final_rwm1[-1,],band = list(number = 8, level = 80))
prctilemlx(final_rwm2[-1,],band = list(number = 8, level = 80))
prctilemlx(final_rwm3[-1,],band = list(number = 8, level = 80))

final_new <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_linear<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))
  post_new<-mcmc(saemix.model,saemix.data,saemix.options_linear)$eta[[indiv]]
  post_new['individual'] <- j
  final_new <- rbind(final_new,post_new)
}




names(final_new)[1]<-paste("time")
names(final_new)[5]<-paste("id")
final_new1 <- final_new[c(5,1,2)]
final_new2 <- final_new[c(5,1,3)]
final_new3 <- final_new[c(5,1,4)]


prctilemlx(final_new1[-1,],band = list(number = 8, level = 80))
prctilemlx(final_new2[-1,],band = list(number = 8, level = 80))
prctilemlx(final_new3[-1,],band = list(number = 8, level = 80))



