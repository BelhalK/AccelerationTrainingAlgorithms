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



K1 = 50
K2 = 20
iterations = 1:(K1+K2+1)
gd_step = 0.00001
replicate = 10
seed0 = 39546

#RWM
final_rwm <- 0
for (j in 1:replicate){
  print(j)
  options<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0)
  theo_ref<-data.frame(saemix(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- j
  final_rwm <- rbind(final_rwm,theo_ref)
}



names(final_rwm)[1]<-paste("time")
names(final_rwm)[9]<-paste("id")
final_rwm1 <- final_rwm[c(9,1,3)]
# prctilemlx(final_rwm1[-1,],band = list(number = 8, level = 80)) + ggtitle("RWM")

#mix (RWM and MAP new kernel for liste of saem iterations)
final_mix <- 0
for (j in 1:replicate){
  print(j)
  options.mix<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,4,0),nbiter.saemix = c(K1,K2),step.gd=gd_step,map.range=c(1:3))
  theo_mix<-data.frame(saemix_gd_mix(saemix.model,saemix.data,options.mix))
  theo_mix <- cbind(iterations, theo_mix)
  theo_mix['individual'] <- j
  final_mix <- rbind(final_mix,theo_mix)
}



names(final_mix)[1]<-paste("time")
names(final_mix)[9]<-paste("id")
final_mix1 <- final_mix[c(9,1,3)]
# prctilemlx(final_mix1[-1,],band = list(number = 8, level = 80)) + ggtitle("mix")

# #map always 
# final_map <- 0
# for (j in 1:replicate){
#   print(j)
#   options.new<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5),nbiter.saemix = c(K1,K2))
#   theo_new_ref<-data.frame(saemix_new(saemix.model,saemix.data,options.new))
#   theo_new_ref <- cbind(iterations, theo_new_ref)
#   theo_new_ref['individual'] <- j
#   final_map <- rbind(final_map,theo_new_ref)
# }


# names(final_map)[1]<-paste("time")
# names(final_map)[9]<-paste("id")
# final_map1 <- final_map[c(9,1,3)]
# prctilemlx(final_map1[-1,],band = list(number = 8, level = 80)) + ggtitle("map")

# for(i in 2:7){
#   final_rwm1 <- final_rwm[c(9,1,i)]
#   final_map1 <- final_map[c(9,1,i)]
  


#  final_rwm1['group'] <- 1
# final_map1['group'] <- 2
# final_map1$id <- final_map1$id +2


#   final <- 0
#   final <- rbind(final_rwm1[-1,],final_map1[-1,])



#   labels <- c("rwm","MAP")
#   final <- final[c(1,4,2,3)]
#   perc <- prctilemlx(final, band = list(number = 2, level = 80),group='group', label = labels) + theme(legend.position = "none")+ ggtitle(colnames(final)[4])
#   ggsave(plot = perc, file = paste("file_map_yield",i,".pdf",sep=""))
# }




final_rwm1['group'] <- 1
final_mix1['group'] <- 2
final_mix1$id <- final_mix1$id +1
# final_map1['group'] <- 2
# final_map1$id <- final_map1$id +2

final <- 0
final <- rbind(final_rwm1[-1,],final_mix1[-1,])
# final <- rbind(final_rwm1[-1,],final_mix1[-1,],final_map1[-1,])



labels <- c("rwm","mix")
# labels <- c("rwm","mix","map")
final <- final[c(1,4,2,3)]
prctilemlx(final, band = list(number = 2, level = 80),group='group', label = labels) + theme(legend.position = "none")

