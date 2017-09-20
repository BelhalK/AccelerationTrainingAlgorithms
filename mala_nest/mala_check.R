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

library("mlxR")
library("psych")
library("coda")
library("Matrix")

require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
iter_mcmc = 800
burn = 400
# Doc
theo.saemix<-read.table( "data/theo.saemix.tab",header=T,na=".")
# l <- c(4.02,4.4,4.53,4.4,5.86,4,4.95,4.53,3.1,5.5,4.92,5.3)
# for (i in 1:12){
#   theo.saemix[(i*10-9):(i*10),'Dose'] = l[i]
# }
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

model1cpt<-function(psi,id,xidep) { 
	dose<-xidep[,1]
	tim<-xidep[,2]  
	ka<-psi[id,1]
	V<-psi[id,2]
	CL<-psi[id,3]
	k<-CL/V
	ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
	return(ypred)
}
# Default model, no covariate
# saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))
saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(10,10,1.05),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))
saemix.options_mala<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc,0,0,0), sigma.val = 0.01, gamma.val = 0.001)
mala <- saemix_mala(saemix.model,saemix.data,saemix.options_mala)
post_mala <- mala$post_mala
rates_mala <- mala$rates


final_mala <- post_mala[[1]]
for (i in 2:length(post_mala)) {
  final_mala <- rbind(final_mala, post_mala[[i]])
}
graphConvMC_new(final_mala, title="MALA")



### plot of average log rejection rate function og log stepsize
sigmas <- seq(from = 0.001, to = 0.01, length.out = 20)
avg <- matrix(NA,ncol=2,nrow=length(sigmas))
colnames(avg) <- c("log stepsize gamma","log rejection rate")
i<-1
for (sigma in sigmas){
  print(i)
  saemix.options_mala<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc,0,0,0), sigma.val = sigma, gamma.val = 0.001)
  mala <- saemix_mala(saemix.model,saemix.data,saemix.options_mala)
  rates_mala <- mala$rates
  avg[i,1] <- log(sigma)
  avg[i,2] <- log(mean(rates_mala[-(1:burn)]))
  i <- i+1
}

coeff <- (avg[2,2]-avg[1,2])/(avg[2,1]-avg[1,1])
plot(avg) + title(coeff)




index = 4
names(post_mala[[index]])[2]<-paste("ka")
names(post_mala[[index]])[3]<-paste("V")
names(post_mala[[index]])[4]<-paste("Cl")


graphConvMC_twokernels(post_mala[[index]],post_mala[[index]], title="RWM vs MALA")



final_mala <- post_mala[[1]]
for (i in 2:length(post_mala)) {
  final_mala <- rbind(final_mala, post_mala[[i]])
}
graphConvMC_new(final_mala, title="MALA")
