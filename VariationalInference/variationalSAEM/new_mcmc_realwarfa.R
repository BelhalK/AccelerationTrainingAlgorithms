load("quantvarinf.RData")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
library(rstan)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
# save.image("quantvarinf.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/VariationalInference/variationalSAEM/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('estep_mcmc.R')
  source('variationalinference.R')
  source('main.R')
  source('indiv_vi.R') 
  source('mcmc_final.R')
  source('main_estep.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  source('graphplot.R')

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/VariationalInference/variationalSAEM")


warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/VariationalInference/variationalSAEM/data/warfarin_data.txt", header=T)
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

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.7,7.51,0.0178),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(0.2,0,0,0,0.18,0,0,0,0.03),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


L_mcmc=10000
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta_ref

options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta


# #Hand made ADVI
# options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=2,nbiter.mcmc = c(0,0,0,1,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# maps<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew)
# Gammamap <- maps$Gamma
# etamap <- maps$map
# K=100
# variational.post.options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.gd = c(K),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),Gamma.laplace=Gamma)
# variational.post<-variational.inference(saemix.model_warfa,saemix.data_warfa,variational.post.options)
# variational.post$mu

# #IMH with maps (new kernel)
# options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=etamap,
#         Gamma = Gammamap,
#         nbiter.mcmc = c(0,0,0,0,6,0),nb.chains=1, nbiter.saemix = c(K1,K2),
#         nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# vi<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfavi)$eta



#RSTAN VB

model <- 'data {
          int<lower=0> N;// Number of observations
          vector[N] time; //predictor
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
          concentration ~ normal(123*beta[1]/(beta[2]*(beta[1]-beta[3]))*(exp(-beta[3]*time)-exp(-beta[1]*time)), pres);
        }'


modelstan <- stan_model(model_name = "warfarin",model_code = model)


#Calculate mu and gamma of ELBO optimization
variational.post.options<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nb.chains=1,
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), modelstan = modelstan)
variational.post<-indiv.variational.inference(saemix.model_warfa,saemix.data_warfa,variational.post.options)
mu.vi <- variational.post$mu
Gamma.vi <- variational.post$Gamma


#Independent sampler with ADVI outputs
i <- 10
# test <- etamap
# test[i,] <- etamap[i,] +0.05
eta.vi <- etamap
Gammavi <- Gammamap
eta.vi[i,] <- mu.vi
Gammavi[[i]] <- Gamma.vi
options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=eta.vi,Gamma = Gammavi,
        nbiter.mcmc = c(0,0,0,0,6,0),nb.chains=1, nbiter.saemix = c(K1,K2),
        nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
vi<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfavi)$eta

# #using the samples from vb (drawn from candidate KL posterior)
# options_warfavb<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,1,0,1),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), modelstan = modelstan)
# vb<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfavb)$eta

#Autocorrelation
rwm.obj <- as.mcmc(ref[[10]])
autocorr.plot(rwm.obj[,1]) + title("RWM SAEM Autocorrelation")

new.obj <- as.mcmc(new[[10]])
autocorr.plot(new.obj[,1]) + title("Laplace SAEM Autocorrelation")

# vb.obj <- as.mcmc(vb[[10]])
# autocorr.plot(vb.obj[,1]) + title("vb SAEM Autocorrelation")

vi.obj <- as.mcmc(vi[[10]])
autocorr.plot(vi.obj[,1]) + title("vb SAEM Autocorrelation")

#MSJD
mssd(ref[[10]][,1])
mssd(new[[10]][,1])
mssd(vb[[10]][,1])
mssd(vi[[10]][,1])



start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 3))

# #one invdiv
# i = 10
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

# # indetabarvb <- vb[[i]]
# # indexpecvb <- data.frame(apply(indetabarvb[-(1:start_interval),], 2, cummean))
# # indexpecvb$iteration <- 1:(L_mcmc-start_interval)


# # indsdvb <- 0
# # indvar <- data.frame(apply(vb[[i]][-(1:start_interval),]^2, 2, cummean))
# # indmeansq <- data.frame(apply(vb[[i]][-(1:start_interval),], 2, cummean))^2
# # indsdvb <- indsdvb + sqrt(pmax(zero,indvar - indmeansq))
# # indsdvb$iteration <- 1:(L_mcmc-start_interval)


# indetabarvi <- vi[[i]]
# indexpecvi <- data.frame(apply(indetabarvi[-(1:start_interval),], 2, cummean))
# indexpecvi$iteration <- 1:(L_mcmc-start_interval)


# indsdvi <- 0
# indvar <- data.frame(apply(vi[[i]][-(1:start_interval),]^2, 2, cummean))
# indmeansq <- data.frame(apply(vi[[i]][-(1:start_interval),], 2, cummean))^2
# indsdvi <- indsdvi + sqrt(pmax(zero,indvar - indmeansq))
# indsdvi$iteration <- 1:(L_mcmc-start_interval)

# # plotmcmc(indexpecref[,c(4,1:3)],indexpecnew[,c(4,1:3)],title=paste("mean",i))
# # plotconv3(indexpecref[,c(4,1:3)],indexpecnew[,c(4,1:3)],indexpecvb[,c(4,1:3)],title="mean")
# plotconv3(indexpecref[,c(4,1:3)],indexpecnew[,c(4,1:3)],indexpecvi[,c(4,1:3)],title="mean")




#quantiles

i <- 10
qref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[[i]][1:k,dim], 0.05)
    qref[[dim]][k,2] <- quantile(ref[[i]][1:k,dim], 0.5)
    qref[[dim]][k,3] <- quantile(ref[[i]][1:k,dim], 0.95)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], 0.05)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], 0.5)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], 0.95)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(4,1:3)],qnew[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
}




# qvb <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
# for (dim in 1:3){
#   print(dim)
#   for (k in 1:L_mcmc){
#     qvb[[dim]][k,1] <- quantile(vb[[i]][1:k,dim], 0.05)
#     qvb[[dim]][k,2] <- quantile(vb[[i]][1:k,dim], 0.5)
#     qvb[[dim]][k,3] <- quantile(vb[[i]][1:k,dim], 0.95)
#   }
#   qvb[[dim]]$iteration <- 1:L_mcmc
#   # plotmcmc(qref[[dim]][,c(4,1:3)],qvb[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
# }



qvi <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qvi[[dim]][k,1] <- quantile(vi[[i]][1:k,dim], 0.05)
    qvi[[dim]][k,2] <- quantile(vi[[i]][1:k,dim], 0.5)
    qvi[[dim]][k,3] <- quantile(vi[[i]][1:k,dim], 0.95)
  }
  qvi[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(4,1:3)],qvb[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=1)+
      xlab("")+scale_x_log10()+ theme_bw() +ylab(names(df[j]))+ theme(axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=15, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=20)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
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


q1new <- data.frame(cbind(iteration,qnew[[1]][,1],qnew[[2]][,1],qnew[[3]][,1]))
q2new <- data.frame(cbind(iteration,qnew[[1]][,2],qnew[[2]][,2],qnew[[3]][,2]))
q3new <- data.frame(cbind(iteration,qnew[[1]][,3],qnew[[2]][,3],qnew[[3]][,3]))
q1new$quantile <- 1
q2new$quantile <- 2
q3new$quantile <- 3
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])

colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
# plotquantile(quantref,quantnew)



# q1vb <- data.frame(cbind(iteration,qvb[[1]][,1],qvb[[2]][,1],qvb[[3]][,1]))
# q2vb <- data.frame(cbind(iteration,qvb[[1]][,2],qvb[[2]][,2],qvb[[3]][,2]))
# q3vb <- data.frame(cbind(iteration,qvb[[1]][,3],qvb[[2]][,3],qvb[[3]][,3]))
# q1vb$quantile <- 1
# q2vb$quantile <- 2
# q3vb$quantile <- 3
# quantvb <- rbind(q1vb[-c(1:burn),],q2vb[-c(1:burn),],q3vb[-c(1:burn),])


# colnames(quantvb)<-c("iteration","ka","V","k","quantile")



q1vi <- data.frame(cbind(iteration,qvi[[1]][,1],qvi[[2]][,1],qvi[[3]][,1]))
q2vi <- data.frame(cbind(iteration,qvi[[1]][,2],qvi[[2]][,2],qvi[[3]][,2]))
q3vi <- data.frame(cbind(iteration,qvi[[1]][,3],qvi[[2]][,3],qvi[[3]][,3]))
q1vi$quantile <- 1
q2vi$quantile <- 2
q3vi$quantile <- 3
quantvarinf <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])
colnames(quantvarinf)<-c("iteration","ka","V","k","quantile")



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
      xlab("")+scale_x_log10()+ theme_bw() +ylab(names(df[j]))+ theme(axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=15, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=20)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}

# plotquantile3(quantref,quantnew,quantvb)
plotquantile3(quantref,quantnew,quantvarinf)

# gelman.plot(mcmc.list(as.mcmc(ref[[10]])), bin.width = 10, max.bins = 50,confidence = 0.95, transform = FALSE, autoburnin=TRUE, auto.layout = TRUE)
# geweke.plot(mcmc.list(as.mcmc(ref[[10]])), frac1=0.1, frac2=0.5)
# geweke.plot(mcmc.list(as.mcmc(new[[10]])), frac1=0.1, frac2=0.5)


q1vi[,2] <- q1vi[,2] + 1
q2vi[,2] <- q2vi[,2] + 1
q3vi[,2] <- q3vi[,2] + 1

q1vi[,3] <- q1vi[,3] + 8 
q2vi[,3] <- q2vi[,3] + 8
q3vi[,3] <- q3vi[,3] + 8

q1vi[,4] <- q1vi[,4] + 0.01 
q2vi[,4] <- q2vi[,4] + 0.01
q3vi[,4] <- q3vi[,4] + 0.01

