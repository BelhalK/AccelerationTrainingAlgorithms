
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library("rCMA")
load("maps_pk.RData")
# save.image("maps_pk.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Desktop/csda_new/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('main.R')
  source('main_maps.R')
  source('main_estep.R')
  source('main_estepmaps.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  # source('SaemixRes_c.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  
setwd("/Users/karimimohammedbelhal/Desktop/csda_new")
source('graphplot.R') 


warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/csda_new/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")
n <- length(unique(warfa_data$id))
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

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,5,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


##RUNS

K1 = 100
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
warfarun<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
warfarun<-saemixmaps(saemix.model_warfa,saemix.data_warfa,options_warfa)
warfa <- warfarun$par
maps <- warfarun$maps
warfa <- cbind(iterations, warfa)


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(2:5))
# warfanew<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew))
warfarunnew<-saemixmaps(saemix.model_warfa,saemix.data_warfa,options_warfanew)
warfanew <- warfarunnew$par
mapsnew <- warfarunnew$maps 
warfanew <- cbind(iterations, warfanew)


for (i in 10:15) {
dmaps<- data.frame(maps[[i]])
dmaps <- cbind(c(1:9,K1),maps[[i]])
dmapsnew<- data.frame(mapsnew[[i]])
dmapsnew <- cbind(c(1:9,K1),mapsnew[[i]])
# plotconv(dmaps,dmapsnew)
plotconv4(dmaps[-1,],dmapsnew[-1,],warfa[c(2:9,K1),1:4],warfanew[c(2:9,K1),1:4])
}


mapmean <- data.frame(1/n*Reduce("+",maps))
mapmean <- cbind(c(1:9,K1),mapmean)
mapmeannew <- data.frame(1/n*Reduce("+",mapsnew))
mapmeannew <- cbind(c(1:9,K1),mapmeannew)

plotconv4(mapmean[-1,],mapsnew[-1,],warfa[c(2:9,K1),1:4],warfanew[c(2:9,K1),1:4])
plotconv(mapmean[-1,],mapmeannew[-1,])



plotconv3(warfa[-1,],warfanew[-1,],warfanew2[-1,])

plotconv(warfa[,c(1,3,6)],warfanew[,c(1,3,6)])
plotconv(warfa[K1:end,c(1,3,6)],warfanew[K1:end,c(1,3,6)])

options_warfanew2<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(K1:(end)))
warfanew2<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew2))
warfanew2 <- cbind(iterations, warfanew2)

plotconv3(warfa[,c(1,3,6)],warfanew[,c(1,3,6)],warfanew2[,c(1,3,6)])
plotconv3(warfa[,c(1,2,5)],warfanew[,c(1,2,5)],warfanew2[,c(1,2,5)])
plotconv3(warfa[,c(1,4,7)],warfanew[,c(1,4,7)],warfanew2[,c(1,4,7)])
plotconv3(warfa[K1:end,c(1,3,6)],warfanew[K1:end,c(1,3,6)],warfanew2[K1:end,c(1,3,6)])
plotconv3(warfa[K1:end,c(1,2,5)],warfanew[K1:end,c(1,2,5)],warfanew2[K1:end,c(1,2,5)])
plotconv3(warfa[K1:end,c(1,4,7)],warfanew[K1:end,c(1,4,7)],warfanew2[K1:end,c(1,4,7)])

plotconv(warfa[,c(1,3,6)],warfanew2[,c(1,3,6)])
plotconv(warfa[(K1-10):end,c(1,3,6)],warfanew2[(K1-10):end,c(1,3,6)])


replicate = 3
seed0 = 395246

#RWM
final_rwm <- 0
final_new <- 0
final_new2 <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(1,5,1,0,0,0),c(0.8,10,0.8,0,0,0),c(1.2,12,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))

  saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
warfa <- cbind(iterations, warfa)
warfa['individual'] <- m
final_rwm <- rbind(final_rwm,warfa)


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(1:50))
warfanew<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew))
warfanew <- cbind(iterations, warfanew)
warfanew['individual'] <- m
final_new <- rbind(final_new,warfanew)

options_warfanew2<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(K1:(K1+100)))
warfanew2<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew2))
warfanew2 <- cbind(iterations, warfanew2)
warfanew2['individual'] <- m
final_new2 <- rbind(final_new2,warfanew2)

}

graphConvMC_new(final_rwm, title="RWM")
graphConvMC_diff(final_rwm[,c(1,3,6,9)],final_new[,c(1,3,6,9)],final_new2[,c(1,3,6,9)], title="RWM")
graphConvMC_diffzoom(final_rwm[,c(1,3,6,9)],final_new[,c(1,3,6,9)],final_new2[,c(1,3,6,9)],(K1-50),(K1+150), title="RWM")


