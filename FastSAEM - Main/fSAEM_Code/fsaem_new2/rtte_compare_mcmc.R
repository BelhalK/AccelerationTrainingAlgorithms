library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
# save.image("rtte_mcmc_conv.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_new2/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 

  source('main.R')
  source('main_estep.R')
  source('estep_mcmc.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  # source('SaemixRes_c.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  
setwd("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_new2")
source('graphplot.R') 
load("rtte_mcmc_conv.RData")
###RTTE
timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_new2/data/rtte_data.csv", header=T, sep=",")
# timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/data/rttellis.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
saemix.data_rtte<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
n <- length(unique(timetoevent.saemix$id))
timetoevent.model<-function(psi,id,xidep) {
T<-xidep[,1]
y<-xidep[,2]
N <- nrow(psi)
Nj <- length(T)
censoringtime = 20
lambda <- psi[id,1]
beta <- psi[id,2]
init <- which(T==0)
cens <- which(T==censoringtime)
ind <- setdiff(1:Nj, append(init,cens))
hazard <- (beta/lambda)*(T/lambda)^(beta-1)
H <- (T/lambda)^beta
logpdf <- rep(0,Nj)
logpdf[cens] <- -H[cens] + H[cens-1]
logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
return(logpdf)
}

saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))


##RUNS

K1 = 200
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

#Weibull
options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
rtte <- cbind(iterations, rtte)


options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5))
rttenew<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenew))
rtte <- cbind(iterations, rtte)



saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(10.17122,4.577724),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(0.3,0,0,0.3),ncol=2, 
  byrow=TRUE))


L_mcmc=5000
options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rtte)$eta_ref

options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rttenew)$eta


start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 2))

etabarref <- 1/n*Reduce("+",ref)
expecref <- data.frame(apply(etabarref[-(1:start_interval),], 2, cummean))
expecref$iteration <- 1:(L_mcmc-start_interval)


sdref <- 0
for (i in 1:n){
  var <- data.frame(apply(ref[[i]][-(1:start_interval),]^2, 2, cummean))
  meansq <- data.frame(apply(ref[[i]][-(1:start_interval),], 2, cummean))^2
  sdref <- sdref + sqrt(pmax(zero,var - meansq))
}


sdref <- 1/n*sdref
sdref$iteration <- 1:(L_mcmc-start_interval)


etabarnew <- 1/n*Reduce("+",new)
expecnew <- data.frame(apply(etabarnew[-(1:start_interval),], 2, cummean))
expecnew$iteration <- 1:(L_mcmc-start_interval)


sdnew <- 0
for (i in 1:n){
  var <- data.frame(apply(new[[i]][-(1:start_interval),]^2, 2, cummean))
  meansq <- data.frame(apply(new[[i]][-(1:start_interval),], 2, cummean))^2
  sdnew <- sdnew + sqrt(pmax(zero,var - meansq))
}

sdnew <- 1/n*sdnew
sdnew$iteration <- 1:(L_mcmc-start_interval)


plotmcmc(expecref[,c(3,1:2)],expecnew[,c(3,1:2)],title="mean")
plotmcmc(sdref[,c(3,1:2)],sdnew[,c(3,1:2)],title="sd")


etaref <- 1/n*Reduce("+",ref)
etaref$iteration <- 1:(L_mcmc)
# plotmcmc(etaref[,c(3,1:2)],etaref[,c(3,1:2)],title="mean")

etanew <- 1/n*Reduce("+",new)
etanew$iteration <- 1:(L_mcmc)

plotmcmc(etaref[,c(3,1:2)],etanew[,c(3,1:2)],title="mean")



# for (i in 1:5){
# ref[[i]]$iteration <- 1:(L_mcmc)
# new[[i]]$iteration <- 1:(L_mcmc)
# plotmcmc(ref[[i]][,c(3,1:2)],new[[i]][,c(3,1:2)],title="mean")
# }

# plotmcmc(ref[[9]][,c(3,1:2)],new[[9]][,c(3,1:2)],title="mean")
# plotmcmc(ref[[5]][,c(3,1:2)],new[[5]][,c(3,1:2)],title="mean")

# #one invdiv


# start_interval <- 200
# zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 2))


# for (i in 1:3){
# indetabarref <- ref[[i]]
# indexpecref <- data.frame(apply(indetabarref[-(1:start_interval),], 2, cummean))
# indexpecref$iteration <- 1:(L_mcmc-start_interval)


# indsdref <- 0
# indvar <- data.frame(apply(ref[[i]][-(1:start_interval),]^2, 2, cummean))
# indmeansq <- data.frame(apply(ref[[i]][-(1:start_interval),], 2, cummean))^2
# indsdref <- indsdref + sqrt(pmax(zero,indvar - indmeansq))
# indsdref$iteration <- 1:(L_mcmc-start_interval)


# indetabarnew <- new[[i]]
# indexpecnew <- data.frame(apply(indetabarnew[-(1:start_interval),], 2, cummean))
# indexpecnew$iteration <- 1:(L_mcmc-start_interval)


# indsdnew <- 0
# indvar <- data.frame(apply(new[[i]][-(1:start_interval),]^2, 2, cummean))
# indmeansq <- data.frame(apply(new[[i]][-(1:start_interval),], 2, cummean))^2
# indsdnew <- indsdnew + sqrt(pmax(zero,indvar - indmeansq))
# indsdnew$iteration <- 1:(L_mcmc-start_interval)


# plotmcmc(indexpecref[,c(3,1:2)],indexpecnew[,c(3,1:2)],title=paste("mean",i))
# plotmcmc(indsdref[-c(1:10),c(3,1:2)],indsdnew[-c(1:10),c(3,1:2)],title=paste("sd",i))
# }



i <- 2
start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 2))
#mean and sd 

indetabarref <- ref[[i]]
indexpecref <- data.frame(apply(indetabarref[-(1:start_interval),], 2, cummean))
indexpecref$iteration <- 1:(L_mcmc-start_interval)


indsdref <- 0
indvar <- data.frame(apply(ref[[i]][-(1:start_interval),]^2, 2, cummean))
indmeansq <- data.frame(apply(ref[[i]][-(1:start_interval),], 2, cummean))^2
indsdref <- indsdref + sqrt(pmax(zero,indvar - indmeansq))
indsdref$iteration <- 1:(L_mcmc-start_interval)


indetabarnew <- new[[i]]
indexpecnew <- data.frame(apply(indetabarnew[-(1:start_interval),], 2, cummean))
indexpecnew$iteration <- 1:(L_mcmc-start_interval)


indsdnew <- 0
indvar <- data.frame(apply(new[[i]][-(1:start_interval),]^2, 2, cummean))
indmeansq <- data.frame(apply(new[[i]][-(1:start_interval),], 2, cummean))^2
indsdnew <- indsdnew + sqrt(pmax(zero,indvar - indmeansq))
indsdnew$iteration <- 1:(L_mcmc-start_interval)


plotmcmc(indexpecref[,c(3,1:2)],indexpecnew[,c(3,1:2)],title=paste("mean",i))
plotmcmc(indsdref[-c(1:10),c(3,1:2)],indsdnew[-c(1:10),c(3,1:2)],title=paste("sd",i))

#quantiles


qref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[[i]][1:k,dim], 0.1)
    qref[[dim]][k,2] <- quantile(ref[[i]][1:k,dim], 0.5)
    qref[[dim]][k,3] <- quantile(ref[[i]][1:k,dim], 0.9)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], 0.1)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], 0.5)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], 0.9)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(3,1:2)],qnew[[dim]][,c(3,1:2)],title=paste("quantiles",i,"dim", dim))
}


plotquantile <- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    if (j<3){
      grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(lambda[i])))
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    }else{
      grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(beta[i])))
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    }
    
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=2, top=title))
}

iteration <- 1:L_mcmc
burn <- 100
q1ref <- data.frame(cbind(iteration,qref[[1]][,1],qref[[2]][,1]))
q2ref <- data.frame(cbind(iteration,qref[[1]][,2],qref[[2]][,2]))
q3ref <- data.frame(cbind(iteration,qref[[1]][,3],qref[[2]][,3]))
q1ref$quantile <- 1
q2ref$quantile <- 2
q3ref$quantile <- 3
quantref <- rbind(q1ref[-c(1:burn),],q2ref[-c(1:burn),],q3ref[-c(1:burn),])


q1new <- data.frame(cbind(iteration,qnew[[1]][,1],qnew[[2]][,1]))
q2new <- data.frame(cbind(iteration,qnew[[1]][,2],qnew[[2]][,2]))
q3new <- data.frame(cbind(iteration,qnew[[1]][,3],qnew[[2]][,3]))
q1new$quantile <- 1
q2new$quantile <- 2
q3new$quantile <- 3
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])




q1ref[,2] <- q1ref[,2] + 10 
q2ref[,2] <- q2ref[,2] + 10
q3ref[,2] <- q3ref[,2] + 10

q1new[,2] <- q1new[,2] + 10
q2new[,2] <- q2new[,2] + 10
q3new[,2] <- q3new[,2] + 10


q1ref[,3] <- q1ref[,3] + 3 
q2ref[,3] <- q2ref[,3] + 3
q3ref[,3] <- q3ref[,3] + 3

q1new[,3] <- q1new[,3] + 3
q2new[,3] <- q2new[,3] + 3
q3new[,3] <- q3new[,3] + 3



colnames(quantref) <- colnames(quantnew)<-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")
save = plotquantile(quantref,quantnew)

save = grid.arrange(save,ncol=1)
ggsave(save, file="newpics/quant_tte.pdf", width = 900, height = 450, units = "mm")
ggsave(save, file="/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_new2/newpics/quantttenotlog.pdf", width = 900, height = 450, units = "mm")

# geweke.plot(mcmc.list(as.mcmc(ref[[10]])), frac1=0.1, frac2=0.5)
# geweke.plot(mcmc.list(as.mcmc(new[[10]])), frac1=0.1, frac2=0.5)


# #cdf
# cdfref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])

# for (dim in 1:2){
#   print(dim)
#   qf1 <- quantile(ref[[i]][,dim], 0.05)
#   qf2 <- quantile(ref[[i]][,dim], 0.5)
#   qf3 <- quantile(ref[[i]][,dim], 0.95)

#   for (k in 1:L_mcmc){
#     cdfref[[dim]][k,1] <- mean(ref[[i]][which(ref[[i]][1:k,dim] < qf1),dim])
#     cdfref[[dim]][k,2] <- mean(ref[[i]][which(ref[[i]][1:k,dim] < qf2),dim])
#     cdfref[[dim]][k,3] <- mean(ref[[i]][which(ref[[i]][1:k,dim] < qf3),dim])
#   }
#   cdfref[[dim]]$iteration <- 1:L_mcmc
# }


# cdfnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
# for (dim in 1:2){
#   print(dim)
#   qf1 <- quantile(new[[i]][,dim], 0.05)
#   qf2 <- quantile(new[[i]][,dim], 0.5)
#   qf3 <- quantile(new[[i]][,dim], 0.95)

#   for (k in 1:L_mcmc){
#     cdfnew[[dim]][k,1] <- mean(new[[i]][which(new[[i]][1:k,dim] < qf1),dim])
#     cdfnew[[dim]][k,2] <- mean(new[[i]][which(new[[i]][1:k,dim] < qf2),dim])
#     cdfnew[[dim]][k,3] <- mean(new[[i]][which(new[[i]][1:k,dim] < qf3),dim])
#   }
#   cdfnew[[dim]]$iteration <- 1:L_mcmc
# }



# iteration <- 1:L_mcmc

# cdf1ref <- data.frame(cbind(iteration,cdfref[[1]][,1],cdfref[[2]][,1]))
# cdf2ref <- data.frame(cbind(iteration,cdfref[[1]][,2],cdfref[[2]][,2]))
# cdf3ref <- data.frame(cbind(iteration,cdfref[[1]][,3],cdfref[[2]][,3]))
# cdf1ref$quantile <- 1
# cdf2ref$quantile <- 2
# cdf3ref$quantile <- 3
# cdfref <- rbind(cdf1ref[-c(1:10),],cdf2ref[-c(1:10),],cdf3ref[-c(1:10),])



# cdf1new <- data.frame(cbind(iteration,cdfnew[[1]][,1],cdfnew[[2]][,1]))
# cdf2new <- data.frame(cbind(iteration,cdfnew[[1]][,2],cdfnew[[2]][,2]))
# cdf3new <- data.frame(cbind(iteration,cdfnew[[1]][,3],cdfnew[[2]][,3]))
# cdf1new$quantile <- 1
# cdf2new$quantile <- 2
# cdf3new$quantile <- 3
# cdfnew <- rbind(cdf1new[-c(1:10),],cdf2new[-c(1:10),],cdf3new[-c(1:10),])


# plotquantile(cdfref,cdfnew, title= "cdf")
