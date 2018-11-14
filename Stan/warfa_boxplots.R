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
# load("RData/boxplots_warfa.RData")
load("RData/newboxplots_warfa.RData")
# save.image("RData/newboxplots_warfa.RData")
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


# saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(1,5,1),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))

# nchains = 100
# L_mcmc= 600
nchains = 10
L_mcmc= 500
indiv.index <- 1


options_warfa<-list(seed=395,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=indiv.index)
  ref<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta


listofrefchains <- list(ref,ref)
# listofrefchains <- 0
for (m in 1:nchains){
  options_warfa<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=indiv.index)
  ref<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta
  # listofrefchains <- listofrefchains + ref
  listofrefchains[[m]] <- ref
}
# averageref <- listofrefchains/nchains

listofnewchains <- list(ref,ref)
# listofnewchains <- 0
for (m in 1:nchains){
  options_warfanew<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=indiv.index)
  new<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta

  # listofnewchains <- listofnewchains + new
  listofnewchains[[m]] <- new
}
# averagenew <- listofnewchains/nchains
#MALA
i <- 10
listofmalachains <- list(ref,ref)
# listofmalachains <- 0
for (m in 1:nchains){
 options.mala<-list(seed=39546*m,map=F,fim=F,ll.is=F, av=0, sigma.val=0.002
  ,gamma.val=0.1,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = i)
mala<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.mala)$eta
# listofmalachains <- listofmalachains + mala
listofmalachains[[m]] <- mala
}
# averagemala <- listofmalachains/nchains


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



#NUTS
listofnutschains <- list(ref,ref)
# listofnutschains <- 0
for (m in 1:nchains){
options.vi<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
nuts<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.vi)$eta
  # listofnutschains <- listofnutschains + nuts
listofnutschains[[m]] <- nuts
}
# groundtruth <- listofnutschains/nchains

groundtruthchains <- list(ref,ref)
for (m in 1:nchains){
  options.vi<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=10000,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
  groundtruth<-mcmc.indiv(saemix.model_warfa,saemix.data_warfa,options.vi)$eta
  groundtruthchains[[m]] <- groundtruth
}


########################PLOT##############################################################

variable <- 1
niter1 <- 5
niter2 <- 20
niter3 <- 500

refbox <- ref[nchains,]
malabox <- mala[nchains,]
newbox <- new[nchains,]
nutsbox <- nuts[nchains,]
truthbox <- groundtruth[nchains,]


for (iter in list(niter1,niter2,niter3)){
  dim <- which(list(niter1,niter2,niter3)==iter)
   for (j in 1:nchains){
    refbox[j,dim] <- listofrefchains[[j]][iter,variable]
    newbox[j,dim] <- listofnewchains[[j]][iter,variable]
    malabox[j,dim] <- listofmalachains[[j]][iter,variable]
    nutsbox[j,dim] <- listofnutschains[[j]][iter,variable]
    truthbox[j,dim] <- groundtruthchains[[j]][dim(groundtruthchains[[j]])[1],variable]
  }
}


refbox$iterations <- 1:nchains
newbox$iterations <- 1:nchains
malabox$iterations <- 1:nchains
nutsbox$iterations <- 1:nchains
truthbox$iterations <- 1:nchains

refbox$algo <- "RWM"
newbox$algo <- "nlme-IMH"
malabox$algo <- "MALA"
nutsbox$algo <- "NUTS"
truthbox$algo <- "Truth"

df <- rbind(refbox,newbox,malabox,nutsbox,truthbox)
colnames(df) <- c(as.character(niter1), as.character(niter2), as.character(niter3),"iterations","algo")
df.m <- melt(df, id.var = c("iterations","algo"))
colnames(df.m) <- c("ID","group","var","value")

# index <- which(df.m$group=="Truth"&df.m$var%in%list(as.character(niter1), as.character(niter2)))
index <- which(df.m$group=="Truth"&df.m$var%in%list(as.character(niter1), as.character(niter2),as.character(niter3)))
test <- df.m[-index,]

# save <- ggplot(data=test, aes(x=factor(group), y=value, fill=factor(var)))+
# geom_boxplot()+
# labs(fill="Iteration")+
# geom_boxplot(data=test[test$group=="Truth",],
#             aes(x=factor(group), y=value,),fill="yellow")+
# labs(fill="Iteration")+
# theme_bw()+
# theme(legend.title = element_text(size=40),
#   legend.key.size = unit(5,"line"),
#       legend.text = element_text(size=40),
#       axis.text=element_text(size=32), 
#       axis.title=element_text(size=40),
#       panel.border = element_rect(colour = "black", fill=NA, size=2),
#       plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))+
# xlab("")+ylab("") + geom_hline(yintercept = quantile(truthbox[,variable],0.25))+
# geom_hline(yintercept = quantile(truthbox[,variable],0.5))+
# geom_hline(yintercept = quantile(truthbox[,variable],0.75))

ggplot(data=test, aes(x=factor(group), y=value, fill=factor(var)))+
geom_boxplot()+
labs(fill="Iteration")+
geom_boxplot(data=test[test$group=="Truth",],
            aes(x=factor(group), y=value,),fill="yellow")+
labs(fill="Iteration")+
theme_bw()+ylab("ka")+
theme(legend.title = element_text(size=40),
  legend.key.size = unit(5,"line"),
      legend.text = element_text(size=40),
      axis.text=element_text(size=32), 
      axis.title=element_text(size=40),
      panel.border = element_rect(colour = "black", fill=NA, size=2),
      plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))+
xlab("")+ylab("ka") + geom_hline(yintercept = quantile(truthbox[,variable],0.25), col="purple", linetype="dashed",size=1.5)+
geom_hline(yintercept = quantile(truthbox[,variable],0.5), col="purple", linetype="dashed",size=1.5)+
geom_hline(yintercept = quantile(truthbox[,variable],0.75), col="purple", linetype="dashed",size=1.5)+
geom_hline(yintercept = quantile(truthbox[,variable],0), col="grey", linetype="dashed",size=1.5)+
geom_hline(yintercept = quantile(truthbox[,variable],1), col="grey", linetype="dashed",size=1.5)


save <- ggplot(data=test, aes(x=factor(group), y=value, fill=factor(var)))+
geom_boxplot()+
labs(fill="Iteration")+
geom_boxplot(data=test[test$group=="Truth",],
            aes(x=factor(group), y=value,),fill="yellow")+
labs(fill="Iteration")+
theme_bw()+ylab("ka")+
theme(legend.title = element_text(size=40),
  legend.key.size = unit(5,"line"),
      legend.text = element_text(size=40),
      axis.text=element_text(size=48), 
      axis.title=element_text(size=56),
      panel.border = element_rect(colour = "black", fill=NA, size=2),
      plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))+
xlab("")+ylab("ka") + geom_hline(yintercept = quantile(truthbox[,variable],0.25), col="purple", linetype="dashed",size=1.5)+
geom_hline(yintercept = quantile(truthbox[,variable],0.5), col="purple", linetype="dashed",size=1.5)+
geom_hline(yintercept = quantile(truthbox[,variable],0.75), col="purple", linetype="dashed",size=1.5)+
geom_hline(yintercept = quantile(truthbox[,variable],0), col="grey", linetype="dashed",size=1.5)+
geom_hline(yintercept = quantile(truthbox[,variable],1), col="grey", linetype="dashed",size=1.5)


# save <- ggplot(data=df.m[-index,]) + 
#     geom_boxplot( aes(x=factor(group), y=value, fill=factor(var)), position=position_dodge(1)) +
#     labs(fill="Iteration")+
#      theme_bw() + theme(legend.title = element_text(size=40),legend.text = element_text(size=40),axis.text=element_text(size=32), 
#                  axis.title=element_text(size=40),
#                    panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))+
#      xlab("Algorithm")+ylab("")

ggsave(save, file="newpics/test1.pdf", width = 900, height = 450, units = "mm")
ggsave(save, file="newpics/boxplots_warfa_V2.pdf", width = 900, height = 450, units = "mm")
ggsave(save, file="newpics/boxplots_warfa_k2.pdf", width = 900, height = 450, units = "mm")

########################PLOT##############################################################



# for (j in 1:nchains){
#   refbox[j,] <- listofrefchains[[j]][niter,]
#   newbox[j,] <- listofnewchains[[j]][niter,]
#   malabox[j,] <- listofmalachains[[j]][niter,]
#   nutsbox[j,] <- listofnutschains[[j]][niter,]
#   truthbox[j,] <- groundtruthchains[[j]][niter,]
# }

refbox$algo <- "RWM"
newbox$algo <- "IMH"
malabox$algo <- "MALA"
nutsbox$algo <- "NUTS"
truthbox$algo <- "Truth"

df <- rbind(refbox,newbox,malabox,nutsbox,truthbox)
colnames(df) <- c("ka", "V", "k","algo")
df.m <- melt(df, id.var = "algo")
algo = c("RWM", "MALA","IMH","NUTS","Truth")
algo <- factor(algo,levels = c("RWM", "MALA","NUTS","IMH","Truth"))
par(cex.axis=5)
ggplot(data = df.m, aes(x=algo, y=value)) + geom_boxplot() + facet_wrap(~variable,ncol = 3)+ theme_bw() 



# niter <- 6
# algo = c("RWM", "MALA","IMH","NUTS", "Truth")
# boxplot(averageref[1:niter,1],averagemala[1:niter,1],averagenew[1:niter,1],groundtruth[1:niter,1],groundtruth[,1], names=algo) 
# boxplot(averageref[1:niter,2],averagemala[1:niter,2],averagenew[1:niter,2],groundtruth[1:niter,2],groundtruth[,2], names=algo) 
# boxplot(averageref[1:niter,3],averagemala[1:niter,3],averagenew[1:niter,3],groundtruth[1:niter,3],groundtruth[,3], names=algo) 

# averagenuts <- groundtruth[1:L_mcmc,]

# averageref$algo <- "RWM"
# averagenew$algo <- "IMH"
# averagemala$algo <- "MALA"
# averagenuts$algo <- "NUTS"
# groundtruth$algo <- "Truth"


# niter <- 6
# df <- rbind(averageref[1:niter,],averagenew[1:niter,],averagemala[1:niter,],averagenuts[1:niter,],tail(groundtruth,1000))
# colnames(df) <- c("ka", "V", "k","algo")
# df.m <- melt(df, id.var = "algo")

# algo <- factor(algo,levels = c("RWM", "MALA","NUTS","IMH", "Truth"))
# par(cex.axis=5)
# ggplot(data = df.m, aes(x=algo, y=value)) + geom_boxplot() + facet_wrap(~variable,ncol = 3)+ theme_bw() 
# # save <- ggplot(data = df.m, aes(x=algo, y=value)) + geom_boxplot() + facet_wrap(~variable,ncol = 3)+ theme_bw() 
# ggsave(save,file="boxplots_warfa.pdf", width = 900, height = 450, units = "mm")
# ggsave(fig,file="boxplot_warfa.pdf", width = 900, height = 450, units = "mm")


par(mfrow=c(1,3))
par(cex.axis=2)
boxplot(averageref[1:niter,1],averagemala[1:niter,1],averagenew[1:niter,1],averagenuts[1:niter,1],groundtruth[,1], names=algo) 
boxplot(averageref[1:niter,2],averagemala[1:niter,2],averagenew[1:niter,2],averagenuts[1:niter,2],groundtruth[,2], names=algo) 
boxplot(averageref[1:niter,3],averagemala[1:niter,3],averagenew[1:niter,3],averagenuts[1:niter,3],groundtruth[,3], names=algo) 



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

