require(ggplot2)
require(gridExtra)
require(reshape2)
# save.image("pkmarc_temptanh1.RData")
# save.image("pkmarc.RData")
# load("pkmarc.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/Annealing - Main/annealing/R")
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

setwd("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/Annealing - Main/annealing")
source('plots.R') 
# data <- read.csv("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/Annealing - Main/annealing/data/joint_datanew.csv", header=T,sep=",")
# data <- data[data$ytype==1,]

data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/AccelerationTrainingAlgorithms/Annealing - Main/annealing/data/joint.txt", header=T)
saemix.data<-saemixData(name.data=data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("dose","time"),name.response=c("y"), name.X="time")

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

saemix.model<-saemixModel(model=model1cpt,description="pk",type="structural"
  ,psi0=matrix(c(1,10,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE), error.model="combined")

##RUNS

K1 = 500
K2 = 100
iterations = 1:(K1+K2)
end = K1+K2



# a <- 0.01
# b <- 0.9
# delay <- 10
# kappa <- iterations/delay + 3*pi/4
# T0 <- 1
# Temp1 <- tanh(iterations/delay) + (T0 - 2*sqrt(2)/(3*pi)*b)*a^(kappa/delay)+b*sin(kappa)/kappa
# plot(Temp1)



a <- 0.5
b <- 0.5
delay <- 2
kappa <- iterations/delay + 3*pi/4
T0 <- 1
Temp2 <- tanh(iterations/delay) + (T0 - 2*sqrt(2)/(3*pi)*b)*a^(kappa/delay)+b*sin(kappa)/kappa
plot(Temp2)


#pkrin
options<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),
  nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,
  nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1)
pk<-data.frame(saemix(saemix.model,saemix.data,options))
pk <- cbind(iterations, pk[-1,])

options.sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),
  nb.chains=1, nbiter.saemix = c(K1,K2),displayProgress=FALSE,
  nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1, alpha.sa=0.95)
pk.sa<-data.frame(saemix(saemix.model,saemix.data,options.sa))
pk.sa <- cbind(iterations, pk.sa[-1,])


optionsnew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),
 nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,
 displayProgress=FALSE,nbiter.burn =0, an=TRUE,coeff=0.5)
pknew<-data.frame(saemix(saemix.model,saemix.data,optionsnew))
pknew <- cbind(iterations, pknew[-1,])

#black blue red
plot3(pk,pk.sa,pknew)

# optionsnew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, an=TRUE,coeff=2.5)
# pknew2<-data.frame(saemix(saemix.model,saemix.data,optionsnew))
# pknew2 <- cbind(iterations, pknew2[-1,])

# optionsnew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, an=TRUE,coeff=1.5)
# pknew3<-data.frame(saemix(saemix.model,saemix.data,optionsnew))
# pknew3 <- cbind(iterations, pknew3[-1,])
plot5(pk,pk.sa,pknew,pknew2,pknew3)

ka <- runif(replicate, min=0.2, max=2)
V <- runif(replicate, min=1, max=10)
k <- runif(replicate, min=0.2, max=2)


replicate = 3

final.ref <- 0
final.sa <- 0
final.an <- 0
# l = list(c(1,10,1,0,0,0),c(2,5,2,0,0,0),c(3,12,3,0,0,0))

for (m in 1:replicate){
  print(m)
  l <- c(ka[m], V[m], k[m])
  
  saemix.model<-saemixModel(model=model1cpt,description="pk",type="structural"
  ,psi0=matrix(l,ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE), error.model="combined")

  options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1,
   nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, 
   map.range=c(0),an=FALSE,coeff=1)
  pk<-data.frame(saemix(saemix.model,saemix.data,options))
  pk <- cbind(iterations, pk[-1,])
  pk['individual'] <- m
  final.ref <- rbind(final.ref,pk)


  options.sa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1,
   nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, 
   map.range=c(0),an=FALSE,coeff=1, alpha.sa=0.95)
  pk.sa<-data.frame(saemix(saemix.model,saemix.data,options.sa))
  pk.sa <- cbind(iterations, pk.sa[-1,])
  pk.sa['individual'] <- m
  final.sa <- rbind(final.sa,pk.sa)

  optionsnew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nb.chains=1,
   nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, an=TRUE,coeff=0.06)
  pknew<-data.frame(saemix(saemix.model,saemix.data,optionsnew))
  pknew <- cbind(iterations, pknew[-1,])
  pknew['individual'] <- m
  final.an <- rbind(final.an,pknew)

}


diff(final.ref,final.sa,final.an)
# diff(pk,pk.sa,pknew)
# diff(final.ref[final.ref$individual==1,],final.sa[final.sa$individual==1,],final.an[final.an$individual==1,])
# diff(final.ref[final.ref$individual==2,],final.sa[final.sa$individual==2,],final.an[final.an$individual==2,])
# diff(final.ref[final.ref$individual==3,],final.sa[final.sa$individual==3,],final.an[final.an$individual==3,])


