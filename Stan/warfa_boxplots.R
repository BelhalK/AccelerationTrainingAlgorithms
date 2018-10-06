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
load("boxplots_warfa.RData")
# save.image("boxplots_warfa.RData")
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
  source('mcmc_indiv.R')
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



warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/Stan/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")


n <- length(unique(warfa_data$id))
model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]

  ypred<-dose*ka/(V*(ka-k))*(exp(-k*time)-exp(-ka*time))
  return(ypred)
}

# saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))


##RUNS

K1 = 400
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

# #Warfarin
# options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
# warfa <- cbind(iterations, warfa)


# options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(K1))
# warfanew<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew))
# warfanew <- cbind(iterations, warfanew)


# graphConvMC_twokernels(warfa,warfanew)


#compareMCMC

# saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(0.7,7.51,0.0178),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(0.5,0,0,0,0.2,0,0,0,0.03),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,8,0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(0.2,0,0,0,0.18,0,0,0,0.03),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

nchains = 50
L_mcmc=100
indiv.index <- 10



# options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=10000,
#   nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
#   nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=indiv.index)
# groundtruth<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta


# listofrefchains <- list(ref,ref)
listofrefchains <- 0
for (m in 1:nchains){
  options_warfa<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=indiv.index)
  ref<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta

  listofrefchains <- listofrefchains + ref
}
averageref <- listofrefchains/nchains


listofnewchains <- 0
for (m in 1:nchains){
  options_warfanew<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=indiv.index)
  new<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta

  listofnewchains <- listofnewchains + new
}
averagenew <- listofnewchains/nchains


# boxplot(groundtruth)
# boxplot(averagenew, add=TRUE,col=c('mistyrose'))
# boxplot(averageref, add=TRUE,col=c('powderblue'))
# # new <- boxplot(averagenew,col=c('green'))
# # truth <- boxplot(groundtruth,col=c('mistyrose'))
# # ref <- boxplot(averageref,col=c('blue'))


# niter <- 10
# algo = c("RWM", "IMH", "Truth")
# boxplot(averageref[1:niter,1],averagenew[1:niter,1],groundtruth[,1], names=algo) 

# abline(h=quantile(groundtruth[,1],0.1),col="red",lty=2)
# abline(h=quantile(groundtruth[,1],0.5),col="red",lty=2)
# abline(h=quantile(groundtruth[,1],0.9),col="red",lty=2)

# abline(h=quantile(averagenew[1:niter,1],0.1),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,1],0.5),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,1],0.9),col="blue",lty=2)

# abline(h=quantile(averageref[1:niter,1],0.1),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,1],0.5),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,1],0.9),col="green",lty=2)


# algo = c("RWM", "IMH", "Truth")
# boxplot(averageref[1:niter,2],averagenew[1:niter,2],groundtruth[,2], names=algo) 


# abline(h=quantile(groundtruth[,2],0.1),col="red",lty=2)
# abline(h=quantile(groundtruth[,2],0.5),col="red",lty=2)
# abline(h=quantile(groundtruth[,2],0.9),col="red",lty=2)

# abline(h=quantile(averagenew[1:niter,2],0.1),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,2],0.5),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,2],0.9),col="blue",lty=2)

# abline(h=quantile(averageref[1:niter,2],0.1),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,2],0.5),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,2],0.9),col="green",lty=2)

# algo = c("RWM", "IMH", "Truth")
# boxplot(averageref[1:niter,3],averagenew[1:niter,3],groundtruth[,3], names=algo) 


#MALA
i <- 10
listofmalachains <- 0
for (m in 1:nchains){
 options.mala<-list(seed=39546*m,map=F,fim=F,ll.is=F, av=0, sigma.val=0.002
  ,gamma.val=0.1,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = i)
mala<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.mala)$eta


  listofmalachains <- listofmalachains + mala
}
averagemala <- listofmalachains/nchains


# niter <- 10
# algo = c("RWM", "MALA","IMH", "Truth")
# boxplot(averageref[1:niter,1],averagemala[1:niter,1],averagenew[1:niter,1],groundtruth[,1], names=algo) 

# abline(h=quantile(groundtruth[,1],0.1),col="red",lty=2)
# abline(h=quantile(groundtruth[,1],0.5),col="red",lty=2)
# abline(h=quantile(groundtruth[,1],0.9),col="red",lty=2)

# abline(h=quantile(averagenew[1:niter,1],0.1),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,1],0.5),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,1],0.9),col="blue",lty=2)

# abline(h=quantile(averageref[1:niter,1],0.1),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,1],0.5),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,1],0.9),col="green",lty=2)

# abline(h=quantile(averagemala[1:niter,1],0.1),col="yellow",lty=2)
# abline(h=quantile(averagemala[1:niter,1],0.5),col="yellow",lty=2)
# abline(h=quantile(averagemala[1:niter,1],0.9),col="yellow",lty=2)

# boxplot(averageref[1:niter,2],averagemala[1:niter,2],averagenew[1:niter,2],groundtruth[,2], names=algo) 

# abline(h=quantile(groundtruth[,2],0.1),col="red",lty=2)
# abline(h=quantile(groundtruth[,2],0.5),col="red",lty=2)
# abline(h=quantile(groundtruth[,2],0.9),col="red",lty=2)

# abline(h=quantile(averagenew[1:niter,2],0.1),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,2],0.5),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,2],0.9),col="blue",lty=2)

# abline(h=quantile(averageref[1:niter,2],0.1),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,2],0.5),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,2],0.9),col="green",lty=2)

# abline(h=quantile(averagemala[1:niter,2],0.1),col="yellow",lty=2)
# abline(h=quantile(averagemala[1:niter,2],0.5),col="yellow",lty=2)
# abline(h=quantile(averagemala[1:niter,2],0.9),col="yellow",lty=2)


# boxplot(averageref[1:niter,3],averagemala[1:niter,3],averagenew[1:niter,3],groundtruth[,3], names=algo) 

# abline(h=quantile(groundtruth[,3],0.1),col="red",lty=2)
# abline(h=quantile(groundtruth[,3],0.5),col="red",lty=2)
# abline(h=quantile(groundtruth[,3],0.9),col="red",lty=2)

# abline(h=quantile(averagenew[1:niter,3],0.1),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,3],0.5),col="blue",lty=2)
# abline(h=quantile(averagenew[1:niter,3],0.9),col="blue",lty=2)

# abline(h=quantile(averageref[1:niter,3],0.1),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,3],0.5),col="green",lty=2)
# abline(h=quantile(averageref[1:niter,3],0.9),col="green",lty=2)

# abline(h=quantile(averagemala[1:niter,3],0.1),col="yellow",lty=2)
# abline(h=quantile(averagemala[1:niter,3],0.5),col="yellow",lty=2)
# abline(h=quantile(averagemala[1:niter,3],0.9),col="yellow",lty=2)



niter <- 6
algo = c("RWM", "MALA","IMH","NUTS", "Truth")
boxplot(averageref[1:niter,1],averagemala[1:niter,1],averagenew[1:niter,1],groundtruth[1:niter,1],groundtruth[,1], names=algo) 
boxplot(averageref[1:niter,2],averagemala[1:niter,2],averagenew[1:niter,2],groundtruth[1:niter,2],groundtruth[,2], names=algo) 
boxplot(averageref[1:niter,3],averagemala[1:niter,3],averagenew[1:niter,3],groundtruth[1:niter,3],groundtruth[,3], names=algo) 

averagenuts <- groundtruth[1:L_mcmc,]

averageref$algo <- "RWM"
averagenew$algo <- "IMH"
averagemala$algo <- "MALA"
averagenuts$algo <- "NUTS"
groundtruth$algo <- "Truth"


niter <- 6
df <- rbind(averageref[1:niter,],averagenew[1:niter,],averagemala[1:niter,],averagenuts[1:niter,],tail(groundtruth,1000))
colnames(df) <- c("ka", "V", "k","algo")
df.m <- melt(df, id.var = "algo")

algo <- factor(algo,levels = c("RWM", "MALA","NUTS","IMH", "Truth"))
par(cex.axis=5)
save <- ggplot(data = df.m, aes(x=algo, y=value)) + geom_boxplot() + facet_wrap(~variable,ncol = 3)+ theme_bw() 
ggsave(save,file="boxplots_warfa.pdf", width = 900, height = 450, units = "mm")
# ggsave(fig,file="boxplot_warfa.pdf", width = 900, height = 450, units = "mm")

# abline(h=quantile(groundtruth[,3],0.9),col="red",lty=2)
# abline(h=quantile(averagenew[1:niter,3],0.9),col="blue",lty=2)
# abline(h=quantile(averageref[1:niter,3],0.9),col="green",lty=2)
# abline(h=quantile(averagemala[1:niter,3],0.9),col="yellow",lty=2)
# abline(h=quantile(groundtruth[1:niter,3],0.9),col="brown",lty=2)

# abline(h=quantile(groundtruth[,3],0.1),col="red",lty=2)
# abline(h=quantile(averagenew[1:niter,3],0.1),col="blue",lty=2)
# abline(h=quantile(averageref[1:niter,3],0.1),col="green",lty=2)
# abline(h=quantile(averagemala[1:niter,3],0.1),col="yellow",lty=2)
# abline(h=quantile(groundtruth[1:niter,3],0.1),col="brown",lty=2)

# abline(h=quantile(groundtruth[,2],0.5),col="red",lty=2)
# abline(h=quantile(averagenew[1:niter,2],0.5),col="blue",lty=2)
# abline(h=quantile(groundtruth[1:niter,2],0.5),col="brown",lty=2)
# abline(h=quantile(averageref[1:niter,3],0.5),col="green",lty=2)
# abline(h=quantile(averagemala[1:niter,3],0.5),col="yellow",lty=2)

# boxplot(groundtruth)
# boxplot(groundtruth[1:niter,])
# boxplot(averagenew, add=TRUE,col=c('mistyrose'))
# boxplot(averageref, add=TRUE,col=c('powderblue'))


#RSTAN VB

model <- 'data {
          int<lower=0> N;// Number of observations
          vector[N] time; //predictor
          real dose; //predictor
          vector[N] concentration;  //response
          
          real beta1_pop;
          real beta2_pop;
          real beta3_pop;
          real<lower=0> omega_beta1;
          real<lower=0> omega_beta2;
          real<lower=0> omega_beta3;
          real<lower=0>  pres;
        }
        parameters {
          vector<lower=0>[3] beta;
        }
        model {
          //Priors
          beta[1] ~ lognormal( beta1_pop , omega_beta1);
          beta[2] ~ lognormal( beta2_pop , omega_beta2);
          beta[3] ~ lognormal( beta3_pop , omega_beta3);

          concentration ~ normal(dose*beta[1]/(beta[2]*(beta[1]-beta[3]))*(exp(-beta[3]*time)-exp(-beta[1]*time)), pres);
        }'

modelstan <- stan_model(model_name = "warfarin",model_code = model)

#NUTS using rstan



#MALA
i <- 10
listofnutschains <- 0
for (m in 1:nchains){
options.vi<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=100000,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
nuts<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.vi)$eta


  listofnutschains <- listofnutschains + nuts
}
groundtruth <- listofnutschains/nchains


i <- 10
options.vi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=100000,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
groundtruth<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.vi)$eta
# vi<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.vi)$eta


#ADVI for VI post outputs
#Calculate mu and gamma of ELBO optimization
variational.post.options<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nb.chains=1,
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),
  modelstan = modelstan, indiv.index = i)

variational.post<-indiv.variational.inference(saemix.model_warfa,saemix.data_warfa,variational.post.options)
mu.vi <- variational.post$mu
Gamma.vi <- variational.post$Gamma
etamap <- variational.post$map
Gammamap <- variational.post$Gammamap

# #using the output of ADVI (drawn from candidate KL posterior)
# test <- etamap
# # test[i,] <- etamap[i,] +0.01
# test[i,] <- mu.vi
eta.vi <- etamap
Gammavi <- Gammamap
# eta.vi[i,] <- mu.vi
# Gammavi[[i]] <- Gamma.vi
options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=eta.vi,Gamma = Gammavi,
        nbiter.mcmc = c(0,0,0,0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),
        nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index = i)
advi<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfavi)$eta

abs(mu.vi - etamap[i,])/abs(etamap[i,])
norm(Gamma.vi - Gammamap[[i]])/norm(Gammamap[[i]])
abs(Gamma.vi - Gammamap[[i]])/abs(Gammamap[[i]])
(Gamma.vi - Gammamap[[i]])/(Gammamap[[i]])

