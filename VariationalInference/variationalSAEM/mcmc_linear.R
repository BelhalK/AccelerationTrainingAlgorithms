load("vb_linear.RData")
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
# save.image("vb_linear.RData")
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

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/VariationalInference/variationalSAEM")


require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

K1 = 400
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

# Doc
oxboys.saemix<-read.table( "data/ox_synth.csv",header=T,na=".",sep=",")
oxboys.saemix_less <- oxboys.saemix[,]
n <- length(unique(oxboys.saemix_less$id))

saemix.data<-saemixData(name.data=oxboys.saemix_less,header=TRUE,
  name.group=c("id"),name.predictors=c("time"),name.response=c("y"),
  units=list(x="yr",y="cm"))

# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

growth.linear<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (2 columns, base and slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  base<-psi[id,1]
  slope<-psi[id,2]
  f<-base+slope*x
  return(f)
}

saemix.model<-saemixModel(model=growth.linear,description="Linear model",type="structural",
  psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2,byrow=TRUE),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE), 
  error.model="constant")


L_mcmc=1000
#RWM mcmc
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model,saemix.data,options_warfa)$eta_ref

#New kernel mcmc
options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model,saemix.data,options_warfanew)$eta

#Gamma Laplace computation
options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=2,nbiter.mcmc = c(0,0,0,6,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
Gamma<-mcmc(saemix.model,saemix.data,options_warfanew)$Gamma


K=10000
i=10


#Variational posterior parameters estimation
variational.post.options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.gd = c(K),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),Gamma.laplace=Gamma)
variational.post<-indiv.variational.inference(saemix.model,saemix.data,variational.post.options)
mus <- variational.post$mu
muss <- transpose(as.data.frame(matrix(unlist(mus), nrow=length(unlist(mus[1])))))
muss$iteration = 1:(K+1)
plotmcmc(muss[,c(3,1:2)],muss[,c(3,1:2)],title=paste("mean VI output",i))

#MCMC with VI proposal
options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=variational.post$mu,
        Gamma = variational.post$Gamma,
        nbiter.mcmc = c(0,0,0,0,0,0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),
        nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
vi<-mcmc(saemix.model,saemix.data,options_warfavi)$eta


#RSTAN VB

model <- 'data {
          int<lower=0> N;// Number of observations
          vector[N] age; //predictor
          vector[N] height;  //response
          
          real beta1_pop;
          real beta2_pop;
          real<lower=0> omega_beta1;
          real<lower=0> omega_beta2;
          real<lower=0>  pres;
        }
        parameters {
          vector[2] beta;
        }
        model {
          //Priors
          beta[1] ~ normal( beta1_pop , omega_beta1);
          beta[2] ~ normal( beta2_pop , omega_beta2);
          height ~ normal(beta[1] + beta[2] * age, 1);
        }'


modelstan <- stan_model(model_name = "oxboys",model_code = model)
options.vb<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,0,1),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), modelstan = modelstan)
vb<-mcmc(saemix.model,saemix.data,options.vb)$eta

# m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
# f <- sampling(modelstan)
# fit_samples = extract(f)
# fit_samples$y[1]


#Autocorrelation
rwm.obj <- as.mcmc(ref[[10]])
autocorr.plot(rwm.obj[,1]) + title("RWM SAEM Autocorrelation")

new.obj <- as.mcmc(new[[10]])
autocorr.plot(new.obj[,1]) + title("Laplace SAEM Autocorrelation")

vb.obj <- as.mcmc(vb[[10]])
autocorr.plot(vb.obj[,1]) + title("VI SAEM Autocorrelation")


#MSJD
mssd(ref[[10]][,1])
mssd(new[[10]][,1])
mssd(vb[[10]][,1])


#Mean plot
start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 2))

i = 10
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

indetabarvb <- vb[[i]]
indexpecvb <- data.frame(apply(indetabarvb[-(1:start_interval),], 2, cummean))
indexpecvb$iteration <- 1:(L_mcmc-start_interval)


indsdvb <- 0
indvar <- data.frame(apply(vb[[i]][-(1:start_interval),]^2, 2, cummean))
indmeansq <- data.frame(apply(vb[[i]][-(1:start_interval),], 2, cummean))^2
indsdvb <- indsdvb + sqrt(pmax(zero,indvar - indmeansq))
indsdvb$iteration <- 1:(L_mcmc-start_interval)


# indetabarvi <- vi[[i]]
# indexpecvi <- data.frame(apply(indetabarvi[-(1:start_interval),], 2, cummean))
# indexpecvi$iteration <- 1:(L_mcmc-start_interval)


# indsdvi <- 0
# indvar <- data.frame(apply(vi[[i]][-(1:start_interval),]^2, 2, cummean))
# indmeansq <- data.frame(apply(vi[[i]][-(1:start_interval),], 2, cummean))^2
# indsdvi <- indsdvi + sqrt(pmax(zero,indvar - indmeansq))
# indsdvi$iteration <- 1:(L_mcmc-start_interval)

plotmcmc(indexpecref[,c(3,1:2)],indexpecnew[,c(3,1:2)],title=paste("mean",i))
plotconv3(indexpecref[,c(3,1:2)],indexpecnew[,c(3,1:2)],indexpecvb[,c(3,1:2)],title="mean")


#Quantiles plot

#quantiles

i <- 10
qref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[[i]][1:k,dim], 0.05)
    qref[[dim]][k,2] <- quantile(ref[[i]][1:k,dim], 0.5)
    qref[[dim]][k,3] <- quantile(ref[[i]][1:k,dim], 0.95)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:32){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], 0.05)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], 0.5)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], 0.95)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(4,1:3)],qnew[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
}



qvb <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:2){
  print(dim)
  for (k in 1:L_mcmc){
    qvb[[dim]][k,1] <- quantile(vb[[i]][1:k,dim], 0.05)
    qvb[[dim]][k,2] <- quantile(vb[[i]][1:k,dim], 0.5)
    qvb[[dim]][k,3] <- quantile(vb[[i]][1:k,dim], 0.95)
  }
  qvb[[dim]]$iteration <- 1:L_mcmc
  # plotmcmc(qref[[dim]][,c(4,1:3)],qnew2[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
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
≤
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


colnames(quantref) <- colnames(quantnew)<-c("iteration","A","B","quantile")
plotquantile(quantref,quantnew)



q1vb <- data.frame(cbind(iteration,qvb[[1]][,1],qvb[[2]][,1]))
q2vb <- data.frame(cbind(iteration,qvb[[1]][,2],qvb[[2]][,2]))
q3vb <- data.frame(cbind(iteration,qvb[[1]][,3],qvb[[2]][,3]))
q1vb$quantile <- 1
q2vb$quantile <- 2
q3vb$quantile <- 3
quantvb <- rbind(q1vb[-c(1:burn),],q2vb[-c(1:burn),],q3vb[-c(1:burn),])


colnames(quantvb)<-c("iteration","A","B","quantile")

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

plotquantile3(quantref,quantnew,quantvb)