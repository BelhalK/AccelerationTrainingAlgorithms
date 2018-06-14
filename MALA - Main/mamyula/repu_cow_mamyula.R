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

# setwd("/home/belhal.karimi/Desktop/Belhal/mamyula")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mamyula")
source('mala_main.R')
source('main_estep_mala.R')
source('main_mamyula.R')
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

# Doc
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
#   id : vector of indices 
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
iterations = 1:(K1+K2+1)
gd_step = 0.01
replicate = 50
seed0 = 39546


final_rwm <- 0
for (j in 1:replicate){
  print("RWM")
  print(j) 
  options<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,0,0), nbiter.saemix = c(K1,K2))
  theo_ref<-data.frame(saemix(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- j 
  final_rwm <- rbind(final_rwm,theo_ref)
}




names(final_rwm)[1]<-paste("time")
names(final_rwm)[9]<-paste("id")
final_rwm1 <- final_rwm[c(9,1,2)]
# prctilemlx(final_rwm1[-1,],band = list(number = 2, level = 80)) + ggtitle("RWM")


#mix (RWM and MAP new kernel for liste of saem iterations)
final_mala <- 0
for (j in 1:replicate){
  print("MALA")
  print(j)
  options.mala<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5,0,0),nbiter.saemix = c(K1,K2),sigma.val = 0.000001,gamma.val=0.00001)
theo_mala<-data.frame(saemix_mamyula(saemix.model,saemix.data,options.mala))
theo_mala <- cbind(iterations, theo_mala)
  theo_mala['individual'] <- j
  final_mala <- rbind(final_mala,theo_mala)
}



names(final_mala)[1]<-paste("time")
names(final_mala)[9]<-paste("id")
final_mala1 <- final_mala[c(9,1,2)]
# prctilemlx(final_mala1[-1,],band = list(number = 2, level = 80)) + ggtitle("mala")

#map always 
final_mamyula <- 0
for (j in 1:replicate){
  print("MAMYULA")
  print(j)
  options.mamyula<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,0,5,0),nbiter.saemix = c(K1,K2),sigma.val = 0.0001,gamma.val=0.001)
  theo_mamyula<-data.frame(saemix_mamyula(saemix.model,saemix.data,options.mamyula))
  theo_mamyula <- cbind(iterations, theo_mamyula)
  theo_mamyula['individual'] <- j
  final_mamyula <- rbind(final_mamyula,theo_mamyula)
}


names(final_mamyula)[1]<-paste("time")
names(final_mamyula)[9]<-paste("id")
final_mamyula1 <- final_mamyula[c(9,1,2)]
# prctilemlx(final_mamyula1[-1,],band = list(number = 8, level = 80)) + ggtitle("mamyula")
for(i in 2:7){
  final_rwm1 <- final_rwm[c(9,1,i)]
  final_mala1 <- final_mala[c(9,1,i)]
  final_mamyula1 <- final_mamyula[c(9,1,i)]


  final_rwm1['group'] <- 1
  final_mala1['group'] <- 2
  final_mala1$id <- final_mala1$id +1
  final_mamyula1['group'] <- 3
  final_mamyula1$id <- final_mamyula1$id +2


  final <- 0
  final <- rbind(final_rwm1[-1,],final_mala1[-1,],final_mamyula1[-1,])



  labels <- c("rwm","mala","mamyula")
  final <- final[c(1,4,2,3)]
  perc <- prctilemlx(final, band = list(number = 2, level = 80),group='group', label = labels) + theme(legend.position = "none")+ ggtitle(colnames(final)[4])
  ggsave(plot = perc, file = paste("file_cow",i,".pdf",sep=""))
}





# final_rwm1['group'] <- 1
# final_mala1['group'] <- 2
# final_mala1$id <- final_mala1$id +1
# final_mamyula1['group'] <- 3
# final_mamyula1$id <- final_mamyula1$id +2

# final <- 0
# final <- rbind(final_rwm1[-1,],final_mala1[-1,],final_mamyula1[-1,])



# labels <- c("rwm","mala","mamyula")
# final <- final[c(1,4,2,3)]
# # prctilemlx(final, band = list(number = 2, level = 80),group='group', label = labels) + theme(legend.position = "none")

# # perc <- prctilemlx(final, band = list(number = 2, level = 80),group='group', label = labels) + theme(legend.position = "none")
# # ggsave('percentile_mamyula.png', perc)




for(i in 2:7){
  final_rwm1 <- final_rwm[c(9,1,i)]
  final_mala1 <- final_mala[c(9,1,i)]
  final_mamyula1 <- final_mamyula[c(9,1,i)]


  final_rwm1['group'] <- 1
  final_mala1['group'] <- 2
  final_mala1$id <- final_mala1$id +1
  final_mamyula1['group'] <- 3
  final_mamyula1$id <- final_mamyula1$id +2


  final <- 0
  final <- rbind(final_rwm1[-1,],final_mala1[-1,],final_mamyula1[-1,])



  labels <- c("rwm","mala","mamyula")
  final <- final[c(1,4,2,3)]
  perc <- prctilemlx(final, band = list(number = 2, level = 80),group='group', label = labels) + theme(legend.position = "none")+ ggtitle(colnames(final)[4])
  ggsave(plot = perc, file = paste("file_cow",i,".pdf",sep=""))
}






