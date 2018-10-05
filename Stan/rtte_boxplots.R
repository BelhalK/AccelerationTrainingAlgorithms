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
# load("boxplot_rtte.RData")
# save.image("boxplot_rtte.RData")
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
  source('mcmc_indiv.R')
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

K1 = 200
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

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

nchains = 3
L_mcmc=100
indiv.index <- 2



# options_rtte<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=20000,nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
#   nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),indiv.index = indiv.index)
# groundtruth<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rtte)$eta


# listofrefchains <- list(ref,ref)
listofrefchains <- 0
for (m in 1:nchains){
  options_rtte<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
    nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),indiv.index = indiv.index)
  ref<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rtte)$eta
  listofrefchains <- listofrefchains + ref[[indiv.index]]
}
averageref <- listofrefchains/nchains


listofnewchains <- 0
for (m in 1:nchains){

  options_rttenew<-list(seed=39546*m,map=F,fim=F,ll.is=F,
    L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1,
     nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, 
     map.range=c(0), indiv.index = indiv.index)
  new<-mcmc(saemix.model_rtte,saemix.data_rtte,options_rttenew)$eta
  listofnewchains <- listofnewchains + new[[indiv.index]]
}
averagenew <- listofnewchains/nchains


# boxplot(groundtruth[[indiv.index]][1000:2000,])
# boxplot(averagenew, add=TRUE,col=c('mistyrose'))
# boxplot(averageref, add=TRUE,col=c('powderblue'))
# # new <- boxplot(averagenew,col=c('green'))
# # truth <- boxplot(groundtruth,col=c('mistyrose'))
# # ref <- boxplot(averageref,col=c('blue'))

# niter <- 200
# algo = c("RWM", "IMH", "Truth")
# boxplot(averageref[1:niter,1],averagenew[1:niter,1],groundtruth[[indiv.index]][,1], names=algo) 

# abline(h=quantile(groundtruth[[indiv.index]][1:niter,1],0.1),col="red",lty=2)
# abline(h=quantile(groundtruth[[indiv.index]][1:niter,1],0.5),col="red",lty=2)
# abline(h=quantile(groundtruth[[indiv.index]][1:niter,1],0.9),col="red",lty=2)

# abline(h=quantile(averagenew[1:niter,1],0.1),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,1],0.5),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,1],0.9),col="blue",lty=2)

# abline(h=quantile(averageref[1:niter,1],0.1),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,1],0.5),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,1],0.9),col="green",lty=2)


# algo = c("RWM", "IMH", "Truth")
# boxplot(averageref[1:niter,2],averagenew[1:niter,2],groundtruth[[indiv.index]][1000:2000,2], names=algo) 


# abline(h=quantile(groundtruth[[indiv.index]][1:niter,2],0.1),col="red",lty=2)
# abline(h=quantile(groundtruth[[indiv.index]][1:niter,2],0.5),col="red",lty=2)
# abline(h=quantile(groundtruth[[indiv.index]][1:niter,2],0.9),col="red",lty=2)

# abline(h=quantile(averagenew[1:niter,2],0.1),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,2],0.5),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,2],0.9),col="blue",lty=2)

# abline(h=quantile(averageref[1:niter,2],0.1),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,2],0.5),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,2],0.9),col="green",lty=2)



# boxplot(averageref, add=TRUE,col=c('powderblue'))
# boxplot(groundtruth, add=TRUE,col=c('mistyrose'))



listofmalachains <- 0
for (m in 1:nchains){
options.mala<-list(seed=39546*m,map=F,fim=F,ll.is=F, av=0, sigma.val=0.01
  ,gamma.val=0.0001,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = indiv.index)
mala<-mcmc(saemix.model_rtte,saemix.data_rtte,options.mala)$eta
  listofmalachains <- listofmalachains + mala[[indiv.index]]
}
averagemala <- listofmalachains/nchains



niter <- 60
algo = c("RWM", "MALA","IMH","NUTS", "Truth")
boxplot(averageref[1:niter,1],averagemala[1:niter,1],averagenew[1:niter,1],averagenuts[1:niter,1],truth[,1], names=algo) 
boxplot(averageref[1:niter,2],averagemala[1:niter,2],averagenew[1:niter,2],groundtruth[1:niter,2],groundtruth[,2], names=algo) 
boxplot(averageref[1:niter,3],averagemala[1:niter,3],averagenew[1:niter,3],groundtruth[1:niter,3],groundtruth[,3], names=algo) 

averagenuts <- data.frame(groundtruth[[indiv.index]])
truth <- data.frame(groundtruth[[indiv.index]])

averageref$algo <- "RWM"
averagenew$algo <- "IMH"
averagemala$algo <- "MALA"
averagenuts$algo <- "NUTS"
truth$algo <- "Truth"

colnames(truth) <-colnames(averageref) <- colnames(averagenew) <- colnames(averagemala) <- colnames(averagenuts) <- c("lambda", "beta","algo")
niter <- 10
# df <- rbind(averageref[1:niter,],averagenew[1:niter,],averagemala[1:niter,],averagemala[1:niter,],averagenuts[1:niter,])
df <- rbind(averageref[1:niter,],averagenew[1:niter,],averagemala[1:niter,],averagenuts[1:niter,], truth)
colnames(df) <- c("lambda", "beta","algo")
df.m <- melt(df, id.var = "algo")

ggplot(data = df.m, aes(x=algo, y=value)) + geom_boxplot() + facet_wrap(~variable,ncol = 3)+ theme_bw() 


# niter <- 1000
# algo = c("IMH", "MALA", "Truth")
# boxplot(averagemala[1:niter,1],averagenew[1:niter,1],groundtruth[[indiv.index]][1000:2000,1], names=algo) 


# algo = c("IMH", "MALA", "Truth")
# boxplot(averagemala[1:niter,2],averagenew[1:niter,2],groundtruth[[indiv.index]][1000:2000,2], names=algo) 




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
options.vi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=10000,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = indiv.index)
groundtruth<-mcmc(saemix.model_rtte,saemix.data_rtte,options.vi)$eta
averagenuts <- groundtruth[[indiv.index]]
truth <- groundtruth[[indiv.index]]
# vi<-mcmc(saemix.model_rtte,saemix.data_rtte,options.vi)$eta



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


