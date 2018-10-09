library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
library(rstan)
# load("rtte_mala_indiv.RData")
load("rtte_mala.RData")
# save.image("rtte_mala.RData")
# save.image("rtte_mala_indiv.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
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

# saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
#   psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
#   c("lambda","beta"))), 
#   transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
#   byrow=TRUE))


##RUNS

# K1 = 200
# K2 = 100
# iterations = 1:(K1+K2+1)
# end = K1+K2

# #Weibull
# options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# rtte<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rtte))
# rtte <- cbind(iterations, rtte)


# options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5))
# rttenew<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options_rttenew))
# rtte <- cbind(iterations, rtte)



saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(10,3),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),omega.init=matrix(c(0.3,0,0,0.3),ncol=2,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

i <- 2
L_mcmc=10000

options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),indiv.index = i)
ref<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rtte)$eta

options_rttenew<-list(seed=39546,map=F,fim=F,ll.is=F,
  L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1,
   nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, 
   map.range=c(0), indiv.index = i)
new<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rttenew)$eta



options.mala<-list(seed=39546,map=F,fim=F,ll.is=F, av=0, sigma.val=0.01
  ,gamma.val=0.0001,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = i)
mala<-mcmc(saemix.model_rtte,saemix.data_rtte,options.mala)$eta


# model <- 'data {
#   int<lower=1> N_e; // Number of total observed events
#   int<lower=1> N_c; // Number of total censoring times ( = # observation units)
#   vector<lower=0>[N_e] event_times; // Times of event occurrence
#   vector<lower=0>[N_c] cens_times; // Censoring times (maybe just N copies of T_c?)
#   real<lower=0> alpha_pop;
#   real<lower=0> sigma_pop;
#   real<lower=0> omega_alpha;
#   real<lower=0> omega_sigma;
# }
# parameters {
#   real<lower=0> alpha; // Weibull shape
#   real<lower=0> sigma; // Weibull scale
# }
# model {
#   // prior
#   alpha ~ lognormal(alpha_pop, omega_alpha);
#   sigma ~ lognormal(sigma_pop, omega_sigma);
  
#   // likelihood
#   for (n_e in 1:N_e) {
#     target += weibull_lpdf(event_times[n_e] | alpha, sigma) - 
#               weibull_lccdf(event_times[n_e] | alpha, sigma);
#   } 
#   for (n_c in 1:N_c) {
#     target += weibull_lccdf(cens_times[n_c] | alpha, sigma);
#   }
# }'

model <- 'data {
  int<lower=1> N_e; // Number of total observed events
  int<lower=1> N_c; // Number of total censoring time
  vector<lower=0>[N_e] event_times; // Times of event occurrence
  int<lower=0> cens_times; // Censoring time
  real<lower=0> beta_pop;
  real<lower=0> lambda_pop;
  real<lower=0> omega_beta;
  real<lower=0> omega_lambda;
}

parameters {
  vector<lower=0>[2] param;
}

model {
  // prior
  param[2] ~ lognormal(beta_pop, omega_beta);
  param[1] ~ lognormal(lambda_pop, omega_lambda);
  
  // likelihood
  target += weibull_lpdf(event_times | param[2], param[1]) - 
            weibull_lccdf(event_times | param[2], param[1]) +
            weibull_lccdf(cens_times | param[2], param[1]);
}'


modelstan <- stan_model(model_name = "rtte",model_code = model)
#NUTS using rstan
options.vi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
vi<-mcmc(saemix.model_rtte,saemix.data_rtte,options.vi)$eta



#ADVI for VI post outputs
#Calculate mu and gamma of ELBO optimization
variational.post.options<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nb.chains=1,
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),
  modelstan = modelstan, indiv.index = i)

variational.post<-indiv.variational.inference(saemix.model_rtte,saemix.data_rtte,variational.post.options)
mu.vi <- variational.post$mu
Gamma.vi <- variational.post$Gamma
etamap <- variational.post$map
Gammamap <- variational.post$Gammamap
# #using the output of ADVI (drawn from candidate KL posterior)

eta.vi <- etamap
Gammavi <- Gammamap
eta.vi[i,] <- mu.vi
Gammavi[[i]] <- Gamma.vi
options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=eta.vi,Gamma = Gammavi,
        nbiter.mcmc = c(0,0,0,0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),
        nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
advi<-mcmc(saemix.model_rtte,saemix.data_rtte,options_warfavi)$eta







par(mfrow=c(1,4))
acf(ref[[i]][,1], main="RWM")
acf(new[[i]][,1], main="IMH")
acf(mala[[i]][,1], main="MALA")
acf(vi[[i]][,1], main="NUTS")

#MSJD
mssd(ref[[i]][,1])
mssd(new[[i]][,1])
mssd(mala[[i]][,1])
mssd(advi[[i]][,1])
mssd(vi[[i]][,1])


i <- 2
start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 2))


#quantiles
qlow <- 0.1
qmed <- 0.5
qhigh <- 0.9


qref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[[i]][1:k,dim], qlow)
    qref[[dim]][k,2] <- quantile(ref[[i]][1:k,dim], qmed)
    qref[[dim]][k,3] <- quantile(ref[[i]][1:k,dim], qhigh)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], qlow)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], qmed)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], qhigh)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(3,1:2)],qnew[[dim]][,c(3,1:2)],title=paste("quantiles",i,"dim", dim))
}

qmala <- list(mala[[i]][1:L_mcmc,],mala[[i]][1:L_mcmc,],mala[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qmala[[dim]][k,1] <- quantile(mala[[i]][1:k,dim], qlow)
    qmala[[dim]][k,2] <- quantile(mala[[i]][1:k,dim], qmed)
    qmala[[dim]][k,3] <- quantile(mala[[i]][1:k,dim], qhigh)
  }
  qmala[[dim]]$iteration <- 1:L_mcmc
}

qvi <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qvi[[dim]][k,1] <- quantile(vi[[i]][1:k,dim], qlow)
    qvi[[dim]][k,2] <- quantile(vi[[i]][1:k,dim], qmed)
    qvi[[dim]][k,3] <- quantile(vi[[i]][1:k,dim], qhigh)
  }
  qvi[[dim]]$iteration <- 1:L_mcmc
}

qadvi <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qadvi[[dim]][k,1] <- quantile(advi[[i]][1:k,dim], qlow)
    qadvi[[dim]][k,2] <- quantile(advi[[i]][1:k,dim], qmed)
    qadvi[[dim]][k,3] <- quantile(advi[[i]][1:k,dim], qhigh)
  }
  qadvi[[dim]]$iteration <- 1:L_mcmc
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


q1mala <- data.frame(cbind(iteration,qmala[[1]][,1],qmala[[2]][,1]))
q2mala <- data.frame(cbind(iteration,qmala[[1]][,2],qmala[[2]][,2]))
q3mala <- data.frame(cbind(iteration,qmala[[1]][,3],qmala[[2]][,3]))
q1mala$quantile <- 1
q2mala$quantile <- 2
q3mala$quantile <- 3
quantmala <- rbind(q1mala[-c(1:burn),],q2mala[-c(1:burn),],q3mala[-c(1:burn),])


q1vi <- data.frame(cbind(iteration,qvi[[1]][,1],qvi[[2]][,1]))
q2vi <- data.frame(cbind(iteration,qvi[[1]][,2],qvi[[2]][,2]))
q3vi <- data.frame(cbind(iteration,qvi[[1]][,3],qvi[[2]][,3]))
q1vi$quantile <- 1
q2vi$quantile <- 2
q3vi$quantile <- 3
quantnuts <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])


q1advi <- data.frame(cbind(iteration,qadvi[[1]][,1],qadvi[[2]][,1]))
q2advi <- data.frame(cbind(iteration,qadvi[[1]][,2],qadvi[[2]][,2]))
q3advi <- data.frame(cbind(iteration,qadvi[[1]][,3],qadvi[[2]][,3]))
q1advi$quantile <- 1
q2advi$quantile <- 2
q3advi$quantile <- 3
quantadvi <- rbind(q1advi[-c(1:burn),],q2advi[-c(1:burn),],q3advi[-c(1:burn),])

colnames(quantref) <- colnames(quantnew)<-colnames(quantmala)<-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")
colnames(quantnuts)<-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")
colnames(quantadvi)<-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")
plotquantile3 <- function(df,df2,df3, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  df3$quantile <- as.factor(df3$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=1)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype = 2,size=1)+
      xlab("")+ theme_bw() +ylab(names(df[j]))+ theme(axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=15, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=20)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=2, top=title))
}

plotquantile3(quantref,quantnew,quantnew)
plotquantile3(quantref,quantnew,quantmala)
plotquantile3(quantref,quantnew,quantnuts)
plotquantile3(quantref,quantnew,quantadvi)