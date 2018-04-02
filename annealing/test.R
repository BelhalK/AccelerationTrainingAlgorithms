
require(ggplot2)
require(gridExtra)
require(reshape2)
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/annealing/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('main.R')
  source('main_estep.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  # source('SaemixRes_c.R') 
  source('SaemixObject.R') 
  source('zzz.R') 

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/annealing")
source('plots.R') 
warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/annealing/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

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
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

##RUNS

K1 = 300
K2 = 100
iterations = 1:(K1+K2)
end = K1+K2

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1)
warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
warfa <- cbind(iterations, warfa[-1,])

options_warfa.sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1)
warfa.sa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa.sa))
warfa.sa <- cbind(iterations, warfa.sa[-1,])


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, an=TRUE,coeff=0.03)
warfanew<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew))
warfanew <- cbind(iterations, warfanew[-1,])

plot3(warfa,warfa.sa,warfanew)



# replicate = 3

# final.ref <- 0
# final.sa <- 0
# final.an <- 0
# for (m in 1:replicate){
#   print(m)
#   l = list(c(1,7,1,0,0,0),c(2,10,2,0,0,0),c(3,12,3,0,0,0))
  
#   saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))

#   options<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1)
#   pk<-data.frame(saemix(saemix.model,saemix.data_warfa,options))
#   pk <- cbind(iterations, pk[-1,])
#   pk['individual'] <- m
#   final.ref <- rbind(final.ref,pk)


#   options.sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1, alpha.sa=0.95)
#   pk.sa<-data.frame(saemix(saemix.model,saemix.data_warfa,options.sa))
#   pk.sa <- cbind(iterations, pk.sa[-1,])
#   pk.sa['individual'] <- m
#   final.sa <- rbind(final.sa,pk.sa)

#   optionsnew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, an=TRUE,coeff=0.03)
#   pknew<-data.frame(saemix(saemix.model,saemix.data_warfa,optionsnew))
#   pknew <- cbind(iterations, pknew[-1,])
#   pknew['individual'] <- m
#   final.an <- rbind(final.an,pknew)

# }

# diff(final.ref,final.sa,final.an)


