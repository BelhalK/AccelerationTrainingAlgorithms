library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
library(data.table)
library(rstan)
# load("hmc_quantile_indiv.RData")
# save.image("../RData/hmc_quantile_indiv_student.RData")
load("../RData/hmc_quantile_indiv_student.RData")
# load("../RData/hmc_quantile_indiv_averagechains.RData")
# save.image("hmc_quantile_indiv.RData")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Stan/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('estep_mcmc.R')
  source('indiv_VI.R')
  source('variationalinferencelinear.R')
  source('main.R')
  source('main_estep.R')
  source('mcmc_final.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('check_linearvslaplace.R')
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  # source('SaemixRes_c.R') 
  source('SaemixObject.R') 
  source('zzz.R')
  source('graphplot.R')

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Stan")

d <- ncol(ref)
i <- 10
start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 3))


#quantiles
qlow <- 0.2
qmed <- 0.5
qhigh <- 0.8


qref <- list(ref[1:L_mcmc,],ref[1:L_mcmc,],ref[1:L_mcmc,])
for (dim in 1:d){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[1:k,dim], qlow)
    qref[[dim]][k,2] <- quantile(ref[1:k,dim], qmed)
    qref[[dim]][k,3] <- quantile(ref[1:k,dim], qhigh)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[1:L_mcmc,],new[1:L_mcmc,],new[1:L_mcmc,])
for (dim in 1:d){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[1:k,dim], qlow)
    qnew[[dim]][k,2] <- quantile(new[1:k,dim], qmed)
    qnew[[dim]][k,3] <- quantile(new[1:k,dim], qhigh)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
}


qnew.student <- list(new.student[1:L_mcmc,],new.student[1:L_mcmc,],new.student[1:L_mcmc,])
for (dim in 1:d){
  print(dim)
  for (k in 1:L_mcmc){
    qnew.student[[dim]][k,1] <- quantile(new.student[1:k,dim], qlow)
    qnew.student[[dim]][k,2] <- quantile(new.student[1:k,dim], qmed)
    qnew.student[[dim]][k,3] <- quantile(new.student[1:k,dim], qhigh)
  }
  qnew.student[[dim]]$iteration <- 1:L_mcmc
}



qvi <- list(vi[1:L_mcmc,],vi[1:L_mcmc,],vi[1:L_mcmc,])
for (dim in 1:d){
  print(dim)
  for (k in 1:L_mcmc){
    qvi[[dim]][k,1] <- quantile(vi[1:k,dim], qlow)
    qvi[[dim]][k,2] <- quantile(vi[1:k,dim], qmed)
    qvi[[dim]][k,3] <- quantile(vi[1:k,dim], qhigh)
  }
  qvi[[dim]]$iteration <- 1:L_mcmc
}


qmala <- list(mala[1:L_mcmc,],mala[1:L_mcmc,],mala[1:L_mcmc,])
for (dim in 1:d){
  print(dim)
  for (k in 1:L_mcmc){
    qmala[[dim]][k,1] <- quantile(mala[1:k,dim], qlow)
    qmala[[dim]][k,2] <- quantile(mala[1:k,dim], qmed)
    qmala[[dim]][k,3] <- quantile(mala[1:k,dim], qhigh)
  }
  qmala[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(4,1:d)],qnew2[[dim]][,c(4,1:d)],title=paste("quantiles",i,"dim", dim))
}

qadvi <- list(advi[1:L_mcmc,],advi[1:L_mcmc,],advi[1:L_mcmc,])
for (dim in 1:d){
  print(dim)
  for (k in 1:L_mcmc){
    qadvi[[dim]][k,1] <- quantile(advi[1:k,dim], qlow)
    qadvi[[dim]][k,2] <- quantile(advi[1:k,dim], qmed)
    qadvi[[dim]][k,3] <- quantile(advi[1:k,dim], qhigh)
  }
  qadvi[[dim]]$iteration <- 1:L_mcmc
}


iteration <- 1:L_mcmc
burn <- 100

q1ref <- data.frame(cbind(iteration,qref[[1]][,1],qref[[2]][,1],qref[[3]][,1]))
q2ref <- data.frame(cbind(iteration,qref[[1]][,2],qref[[2]][,2],qref[[3]][,2]))
q3ref <- data.frame(cbind(iteration,qref[[1]][,3],qref[[2]][,3],qref[[3]][,3]))
q1ref$quantile <- 1
q2ref$quantile <- 2
q3ref$quantile <- 3
quantref <- rbind(q1ref[-c(1:burn),],q2ref[-c(1:burn),],q3ref[-c(1:burn),])
colnames(quantref) <- c("iteration","ka","V","k","quantile")

q1new <- data.frame(cbind(iteration,qnew[[1]][,1],qnew[[2]][,1],qnew[[3]][,1]))
q2new <- data.frame(cbind(iteration,qnew[[1]][,2],qnew[[2]][,2],qnew[[3]][,2]))
q3new <- data.frame(cbind(iteration,qnew[[1]][,3],qnew[[2]][,3],qnew[[3]][,3]))
q1new$quantile <- 1
q2new$quantile <- 2
q3new$quantile <- 3
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])
colnames(quantnew)<-c("iteration","ka","V","k","quantile")


q1new.student <- data.frame(cbind(iteration,qnew.student[[1]][,1],qnew.student[[2]][,1],qnew.student[[3]][,1]))
q2new.student <- data.frame(cbind(iteration,qnew.student[[1]][,2],qnew.student[[2]][,2],qnew.student[[3]][,2]))
q3new.student <- data.frame(cbind(iteration,qnew.student[[1]][,3],qnew.student[[2]][,3],qnew.student[[3]][,3]))
q1new.student$quantile <- 1
q2new.student$quantile <- 2
q3new.student$quantile <- 3
quantnew.student <- rbind(q1new.student[-c(1:burn),],q2new.student[-c(1:burn),],q3new.student[-c(1:burn),])
colnames(quantnew.student)<-c("iteration","ka","V","k","quantile")


q1vi <- data.frame(cbind(iteration,qvi[[1]][,1],qvi[[2]][,1],qvi[[3]][,1]))
q2vi <- data.frame(cbind(iteration,qvi[[1]][,2],qvi[[2]][,2],qvi[[3]][,2]))
q3vi <- data.frame(cbind(iteration,qvi[[1]][,3],qvi[[2]][,3],qvi[[3]][,3]))
q1vi$quantile <- 1
q2vi$quantile <- 2
q3vi$quantile <- 3
quantnuts <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])
colnames(quantnuts)<-c("iteration","ka","V","k","quantile")


q1mala <- data.frame(cbind(iteration,qmala[[1]][,1],qmala[[2]][,1],qmala[[3]][,1]))
q2mala <- data.frame(cbind(iteration,qmala[[1]][,2],qmala[[2]][,2],qmala[[3]][,2]))
q3mala <- data.frame(cbind(iteration,qmala[[1]][,3],qmala[[2]][,3],qmala[[3]][,3]))
q1mala$quantile <- 1
q2mala$quantile <- 2
q3mala$quantile <- 3
quantmala <- rbind(q1mala[-c(1:burn),],q2mala[-c(1:burn),],q3mala[-c(1:burn),])
colnames(quantmala)<-c("iteration","ka","V","k","quantile")


### both VI output
q1advi.full <- data.frame(cbind(iteration,qadvi[[1]][,1],qadvi[[2]][,1],qadvi[[3]][,1]))
q2advi.full <- data.frame(cbind(iteration,qadvi[[1]][,2],qadvi[[2]][,2],qadvi[[3]][,2]))
q3advi.full <- data.frame(cbind(iteration,qadvi[[1]][,3],qadvi[[2]][,3],qadvi[[3]][,3]))
q1advi.full$quantile <- 1
q2advi.full$quantile <- 2
q3advi.full$quantile <- 3
quantadvi.full <- rbind(q1advi.full[-c(1:burn),],q2advi.full[-c(1:burn),],q3advi.full[-c(1:burn),])
colnames(quantadvi.full)<-c("iteration","ka","V","k","quantile")


plotquantile3(quantref,quantnew,quantnew.student)

plotquantile(quantref,quantnew)
plotquantile3(quantref,quantnew,quantmala)
plotquantile3(quantref,quantnew,quantadvi.full)

plotquantile4(quantnew,quantnuts,quantmala,quantadvi.full)

plotquantile4(quantnew,quantnuts,quantnew.student,quantmala)


colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")

par(mfrow=c(3,3))
save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantnew[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantnew[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantnew[,c(1,4,5)])
save <- grid.arrange(save1,save2,save3, ncol=3)

#Autocorrelation
par(mfrow=c(1,3))
acf(ref[,1], main="RWM")
acf(new[,1], main="IMH (Gaussian)")
acf(new.student[,1], main="IMH (Student)")
par(mfrow=c(1,3))
acf(mala[,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(advi[,1], main="ADVI")

par(mfrow=c(1,3))
acf(ref[,2], main="RWM")
acf(new[,2], main="IMH (Gaussian)")
acf(new.student[,2], main="IMH (Student)")
par(mfrow=c(1,3))
acf(mala[,2], main="MALA")
acf(vi[,2], main="NUTS")
acf(advi[,2], main="ADVI")


par(mfrow=c(1,3))
acf(ref[,3], main="RWM")
acf(new[,3], main="IMH (Gaussian)")
acf(new.student[,3], main="IMH (Student)")
par(mfrow=c(1,3))
acf(mala[,3], main="MALA")
acf(vi[,3], main="NUTS")
acf(advi[,3], main="ADVI")




#Autocorrelation
par(mfrow=c(1,6))
acf(ref[,1], main="RWM")
acf(new[,1], main="IMH (Gaussian)")
acf(new.student[,1], main="IMH (Student)")
acf(mala[,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(advi[,1], main="ADVI")

par(mfrow=c(1,6))
acf(ref[,2], main="RWM")
acf(new[,2], main="IMH (Gaussian)")
acf(new.student[,2], main="IMH (Student)")
acf(mala[,2], main="MALA")
acf(vi[,2], main="NUTS")
acf(advi[,2], main="ADVI")


par(mfrow=c(1,6))
acf(ref[,3], main="RWM")
acf(new[,3], main="IMH (Gaussian)")
acf(new.student[,3], main="IMH (Student)")
acf(mala[,3], main="MALA")
acf(vi[,3], main="NUTS")
acf(advi[,3], main="ADVI")






abs(mu.vi - etamap[i,])/abs(etamap[i,])
norm(Gamma.vi - Gammamap[[i]])/norm(Gammamap[[i]])
abs(Gamma.vi - Gammamap[[i]])/abs(Gammamap[[i]])
(Gamma.vi - Gammamap[[i]])/(Gammamap[[i]])

