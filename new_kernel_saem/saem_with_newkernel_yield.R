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

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")


# Doc
data(yield.saemix)
saemix.data<-saemixData(name.data=yield.saemix,header=TRUE,name.group=c("site"),
  name.predictors=c("dose"),name.response=c("yield"),
  name.covariates=c("soil.nitrogen"),units=list(x="kg/ha",y="t/ha",
  covariates=c("kg/ha")))

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



K1 = 100
K2 = 50
iteration = 1:(K1+K2+1)
gd_step = 0.00001


#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0)
theo_ref<-data.frame(saemix(saemix.model,saemix.data,options))
theo_ref <- cbind(iteration, theo_ref)

#ref (map always)
options.new<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5),nbiter.saemix = c(K1,K2))
theo_new_ref<-data.frame(saemix_new(saemix.model,saemix.data,options.new))
theo_new_ref <- cbind(iteration, theo_new_ref)

# #MAP once and  NO GD
# options.nogd<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5),nbiter.saemix = c(K1,K2),step.gd = 0)
# theo_nogd<-data.frame(saemix_gd(saemix.model,saemix.data,options.nogd))
# theo_nogd <- cbind(iterations, theo_nogd)


#MAP once and GD
options.gd<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5),nbiter.saemix = c(K1,K2),step.gd=gd_step)
theo_gd<-data.frame(saemix_gd(saemix.model,saemix.data,options.gd))
theo_gd <- cbind(iteration, theo_gd)

#mix (MAP first 4 iter then gd for 30 and then RWM)
options.gd_mix<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,4),nbiter.saemix = c(K1,K2),step.gd=gd_step)
theo_gd_mix<-data.frame(saemix_gd_mix(saemix.model,saemix.data,options.gd_mix))
theo_gd_mix <- cbind(iteration, theo_gd_mix)

#mix (RWM and MAP new kernel for liste of saem iterations)
options.mix<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,4),nbiter.saemix = c(K1,K2),step.gd=gd_step)
theo_mix<-data.frame(saemix_gd_mix(saemix.model,saemix.data,options.mix))
theo_mix <- cbind(iteration, theo_mix)


theo_mix2<-data.frame(saemix_gd_mix(saemix.model,saemix.data,options.mix))
theo_mix2 <- cbind(iteration, theo_mix2)

graphConvMC_threekernels(theo_ref,theo_mix,theo_mix2, title="ref vs GD")

theo_ref$algo <- 'rwm'
theo_new_ref$algo <- 'MAP'
theo_mix$algo <- 'Mix'

comparison <- 0
comparison <- rbind(theo_ref,theo_mix)
comparison <- rbind(theo_ref,theo_new_ref,theo_mix)
var <- melt(comparison, id.var = c('iteration','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=FALSE)




#RWM vs mix
graphConvMC_twokernels(theo_ref,theo_ref, title="ref vs GD")
graphConvMC_twokernels(theo_ref,theo_new_ref, title="ref vs GD")
graphConvMC_twokernels(theo_ref,theo_gd, title="ref vs GD")
graphConvMC_twokernels(theo_ref,theo_mix, title="ref vs GD")
graphConvMC_threekernels(theo_ref,theo_new_ref,theo_mix, title="ref vs GD")



graphConvMC_twokernels(theo_mix,theo_mix2, title="ref vs GD")
graphConvMC_twokernels(theo_new_ref,theo_mix, title="ref vs GD")



graphConvMC_twokernels(theo_ref,theo_gd_mix, title="ref vs GD")
graphConvMC_twokernels(theo_gd_mix,theo_new_ref, title="ref vs GD")



#RWM vs always MAP (ref)
graphConvMC_twokernels(theo_ref,theo_new_ref, title="new kernel")
#ref vs map once no gd
graphConvMC_twokernels(theo_new_ref,theo_nogd, title="ref vs NOGD")
#map once no gd vs map once and gd
graphConvMC_twokernels(theo_nogd,theo_gd, title="NO GD vs GD")
#ref vs map once gd
graphConvMC_twokernels(theo_new_ref,theo_gd, title="ref vs GD")
#RWM vs GD
graphConvMC_twokernels(theo_ref,theo_gd, title="ref vs GD")

