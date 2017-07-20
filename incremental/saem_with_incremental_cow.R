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
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/incremental")
source('main_incremental.R')
source('main_estep_incremental.R')
source("mixtureFunctions.R")

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")



# Doc
data(cow.saemix)
saemix.data<-saemixData(name.data=cow.saemix,header=TRUE,name.group=c("cow"), 
  name.predictors=c("time"),name.response=c("weight"), 
  name.covariates=c("birthyear","twin","birthrank"), 
  units=list(x="days",y="kg",covariates=c("yr","-","-")))


# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
growthcow<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, a, b, k)
#   id : vector of indices add --all

#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  a<-psi[id,1]
  b<-psi[id,2]
  k<-psi[id,3]
  f<-a*(1-b*exp(-k*x))
  return(f)
}
saemix.model<-saemixModel(model=growthcow,
  description="Exponential growth model", 
  psi0=matrix(c(700,0.9,0.02,0,0,0),ncol=3,byrow=TRUE, 
  dimnames=list(NULL,c("A","B","k"))),transform.par=c(1,1,1),fixed.estim=c(1,1,1), 
  covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
  omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")


K1 = 100
K2 = 50
iteration = 1:(K1+K2+1)



#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE)
theo_ref<-data.frame(saemix(saemix.model,saemix.data,options))
theo_ref <- cbind(iteration, theo_ref)



#ref (map always)
options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),nb.replacement=50,displayProgress=FALSE)
theo_incremental<-data.frame(incremental_saemix(saemix.model,saemix.data,options.incremental))
theo_incremental <- cbind(iteration, theo_incremental)



theo_ref$algo <- 'rwm'
theo_incremental$algo <- 'ISAEM'

theo_ref_scaled <- theo_ref[rep(seq_len(nrow(theo_ref)), each=2),]
theo_ref_scaled$iteration = 1:(2*(K1+K2+1))


comparison <- 0
# comparison <- rbind(theo_ref,theo_incremental)
comparison <- rbind(theo_ref_scaled[iteration,],theo_incremental)

var <- melt(comparison, id.var = c('iteration','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)




