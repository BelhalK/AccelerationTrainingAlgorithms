
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
library(rstan)
require(reshape2)
library("rCMA")

# load("RData/warfa_hmcsaem.RData")
# save.image("RData/warfa_hmcsaem.RData")

# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_newHMC/R")
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
  
setwd("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_newHMC")
source('plotfunc/graphplot.R') 


warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_newHMC/data/warfarin_data.txt", header=T)
saemix.data<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
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

# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE))


# saemix.model<-saemixModel(model=model1cpt,description="warfarin"
#   ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),fixed.estim=c(0,1,0),
#   transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
#   byrow=TRUE))

K1 = 100
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2



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


replicate = 1
seed0 = 395246
seed1 = 3952


final_rwm <- 0
final_mix <- 0

for (m in 1:replicate){
  print(m)
  l = list(c(1,5,1,0,0,0),c(0.8,12,0.8,0,0,0),c(1.2,2,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  # l = list(c(1,5,2,0,0,0),c(3,12,5,0,0,0),c(6,3,7,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),fixed.estim=c(1,1,1),
  transform.par=c(1,1,1),omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE))

  options<-list(seed=seed0/m,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0)
  theo_ref<-data.frame(saemix(saemix.model,saemix.data,options)$par)
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref[,5:7] <- sqrt(theo_ref[,5:7])
  theo_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,theo_ref[-1,])

  options.new<-list(seed=m*seed1,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,2,0),map.range=c(1:2),nbiter.sa=0,nbiter.saemix = c(K1,K2),nbiter.burn =0, modelstan = NULL)
  theo_new_ref<-data.frame(saemix(saemix.model,saemix.data,options.new)$par)
  theo_mix <- cbind(iterations, theo_new_ref)
  theo_mix[,5:7] <- sqrt(theo_mix[,5:7])
  theo_mix['individual'] <- m
  final_mix <- rbind(final_mix,theo_mix[-1,])
}


final_hmc <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(1,5,1,0,0,0),c(0.8,12,0.8,0,0,0),c(1.2,2,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  # l = list(c(1,5,2,0,0,0),c(3,12,5,0,0,0),c(6,3,7,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),fixed.estim=c(1,1,1),
  transform.par=c(1,1,1),omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE))

  options.hmc<-list(seed=seed0/m,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,2),map.range=c(1:2), modelstan = modelstan,nbiter.sa=0,nbiter.burn =0,nbiter.saemix = c(K1,K2))
  theo_hmc_ref<-data.frame(saemix(saemix.model,saemix.data,options.hmc)$par)
  theo_hmc <- cbind(iterations, theo_hmc_ref)
  theo_hmc[,5:7] <- sqrt(theo_hmc[,5:7])
  theo_hmc['individual'] <- m
  final_hmc <- rbind(final_hmc,theo_hmc[-1,])
}

convpop <- graphConvMC_diffvi2(final_rwm[,c(1,3,9)],final_mix[,c(1,3,9)],final_mix[,c(1,3,9)])
convvar <- graphConvMC_diffvi1(final_rwm[,c(1,6,9)],final_mix[,c(1,6,9)],final_mix[,c(1,6,9)])
save <- grid.arrange(convpop,convvar, ncol=2)


convpop <- graphConvMC_diffvi2(final_rwm[,c(1,3,9)],final_mix[,c(1,3,9)],final_hmc[,c(1,3,9)])
convvar <- graphConvMC_diffvi1(final_rwm[,c(1,6,9)],final_mix[,c(1,6,9)],final_hmc[,c(1,6,9)])
save <- grid.arrange(convpop,convvar, ncol=2)


convpop <- graphConvMC_diffvi2(final_rwm[,c(1,2,9)],final_mix[,c(1,2,9)],final_hmc[,c(1,2,9)])
convvar <- graphConvMC_diffvi1(final_rwm[,c(1,5,9)],final_mix[,c(1,5,9)],final_hmc[,c(1,5,9)])
save <- grid.arrange(convpop,convvar, ncol=2)

convpop <- graphConvMC_diffvi2(final_rwm[,c(1,4,9)],final_mix[,c(1,4,9)],final_hmc[,c(1,4,9)])
convvar <- graphConvMC_diffvi1(final_rwm[,c(1,7,9)],final_mix[,c(1,7,9)],final_hmc[,c(1,7,9)])
save <- grid.arrange(convpop,convvar, ncol=2)

# ka_true <- 1
# V_true <- 8
# k_true <- 0.01
# o_ka <- 0.5
# o_V <- 0.2
# o_k <- 0.2
# a_true <- 1

# var_rwm <- 0
# error_rwm <- 0


# var_mix <- 0
# error_mix <- 0

# replicate = 4

# final_rwm <- 0
# final_mix <- 0

# replicate <- 40
# for (m in 1:replicate){
  
  
# myModel <- inlineModel("


# [INDIVIDUAL]
# input = {ka_pop, V_pop, k_pop, omega_ka, omega_V, omega_k}
# DEFINITION:
# ka = {distribution=lognormal, reference=ka_pop, sd=omega_ka}
# V  = {distribution=lognormal, reference=V_pop,  sd=omega_V }
# k = {distribution=lognormal, reference=k_pop, sd=omega_k}


# [LONGITUDINAL]
# input = {ka, V, k,a}
# EQUATION:
# C = pkmodel(ka,V,k)
# DEFINITION:
# y = {distribution=normal, prediction=C, sd=a}
# ")

# N=50

# pop.param   <- c(
#   ka_pop  = ka_true,    omega_ka  = o_ka,
#   V_pop   = V_true,   omega_V   = o_V,
#   k_pop  = k_true,    omega_k  = o_k, a =a_true)
  

# res <- simulx(model     = myModel,
#               parameter = pop.param,
#               treatment = list(time=0, amount=100),
#               group     = list(size=N, level='individual'),
#               output    = list(name='y', time=seq(0,10,by=1)))
  

#   warfarin.saemix <- res$y
#   warfarin.saemix$amount <- 100
  
#   saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(1,10,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE))


# saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("amount","time"),name.response=c("y"), name.X="time")


#   options<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0)
#   theo_ref<-data.frame(saemix(saemix.model,saemix.data,options)$par)
#   theo_ref <- cbind(iterations, theo_ref)
#   theo_ref[,5:8] <- sqrt(theo_ref[,5:8])
#   ML <- theo_ref[,2:8]
#   ML[1:(end+1),]<- theo_ref[end+1,2:8]
#   error_rwm <- error_rwm + (theo_ref[,2:8]-ML)^2
#   theo_ref['individual'] <- m
#   final_rwm <- rbind(final_rwm,theo_ref)
  

#   options.new<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(1:6),nbiter.burn =0)
#   theo_mix<-data.frame(saemix(saemix.model,saemix.data,options.new)$par)
#   theo_mix <- cbind(iterations, theo_mix)
#   theo_mix[,5:8] <- sqrt(theo_mix[,5:8])
#   ML <- theo_mix[,2:8]
#   ML[1:(end+1),]<- theo_mix[end+1,2:8]
#   error_mix <- error_mix + (theo_mix[,2:8]-ML)^2
#   theo_mix['individual'] <- m
#   final_mix <- rbind(final_mix,theo_mix)
# }

# error_rwm <- 1/replicate*error_rwm
# error_mix <- 1/replicate*error_mix

# error_rwm <- cbind(iterations, error_rwm)
# error_mix <- cbind(iterations, error_mix)
# error_rwm[2,] = error_mix[2,]


# err_mix<- theo_ref
# err_rwm<- theo_ref
# err_rwm[,2:8] <- error_rwm[,2:8]
# err_mix[,2:8] <- error_mix[,2:8]


# a <- err_rwm[-1,]
# b <- err_mix[-1,]


# graphConvMC_se1 <- function(df,df2, title=NULL, ylim=NULL)
# {
#   G <- (ncol(df)-2)/3
#   df$individual <- as.factor(df$individual)
#   df2$individual <- as.factor(df2$individual)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
#       xlab("iteration") +scale_x_log10(breaks= c(10,100,200))+ ylab(expression(paste(E(V[pop]))))  
#       grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
#                  axis.title=element_text(size=40),
#                    panel.border = element_rect(colour = "black", fill=NA, size=2))
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj

#   }
#   do.call("grid.arrange", c(graf, ncol=1, top=title))
# }


# graphConvMC_se2 <- function(df,df2, title=NULL, ylim=NULL)
# {
#   G <- (ncol(df)-2)/3
#   df$individual <- as.factor(df$individual)
#   df2$individual <- as.factor(df2$individual)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
#       xlab("iteration") +scale_x_log10(breaks= c(10,100,200))+ ylab(expression(paste(E(omega[V]))))  
#       grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
#                  axis.title=element_text(size=40),
#                    panel.border = element_rect(colour = "black", fill=NA, size=2))
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj

#   }
#   do.call("grid.arrange", c(graf, ncol=1, top=title))
# }


# c <- graphConvMC_se1(a[,c(1,3,9)],b[,c(1,3,9)])
# d <- graphConvMC_se2(a[,c(1,6,9)],b[,c(1,6,9)])


# save <- grid.arrange(c,d, ncol=2)
# # ggsave(save,file="pics_square/sepk.pdf", width = 900, height = 225, units = "mm")

# ggsave(save,file="newpics/se_pk.pdf", width = 900, height = 450, units = "mm")

# # c <- graphConvMC_se1(a[,c(1,2,9)],b[,c(1,2,9)])
# # d <- graphConvMC_se2(a[,c(1,5,9)],b[,c(1,5,9)])

# # grid.arrange(c,d, ncol=2)

# # c <- graphConvMC_se1(a[,c(1,4,9)],b[,c(1,4,9)])
# # d <- graphConvMC_se2(a[,c(1,7,9)],b[,c(1,7,9)])

# # grid.arrange(c,d, ncol=2)


# # graphConvMC_diff4(final_rwm[,c(1,3,9)],final_mix[,c(1,3,9)])