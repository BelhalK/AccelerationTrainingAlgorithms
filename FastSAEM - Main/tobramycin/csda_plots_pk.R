
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library("rCMA")
library(rstan)

  source('R/aaa_generics.R') 
  source('R/compute_LL.R') 
  source('R/func_aux.R') 
  source('R/func_distcond.R') 
  source('R/func_FIM.R')
  source('R/func_plots.R') 
  source('R/func_simulations.R') 
  source('R/main.R')
  source('R/main_estep.R')
  source('R/main_initialiseMainAlgo.R') 
  source('R/main_mstep.R') 
  source('R/SaemixData.R')
  source('R/SaemixModel.R') 
  source('R/SaemixRes.R') 
  # source('R/SaemixRes_c.R') 
  source('R/SaemixObject.R') 
  source('R/zzz.R') 
source('R/graphplot.R') 


# bolus_data <- read.table("data/bolus1_data.txt", header=T)
# saemix.data<-saemixData(name.data=tobra_data,header=TRUE,sep=" ",na=NA, name.group=c("ID"),
#   name.predictors=c("amt","time"),name.response=c("y"), name.X="time")

tobra_data <- read.table("data/tobramycin.txt", header=T)
saemix.data<-saemixData(name.data=tobra_data,header=TRUE,sep=" ",na=NA, name.group=c("ID"),
  name.predictors=c("DOSE","TIME", "EVID"),name.response=c("CP"), name.X="TIME")


model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]
  event<-xidep[,3]  

  dose.index <- which(event==1)
  obs.index <- which(event==0)

  t.dose <- time[dose.index]
  t.obs <- time[obs.index]


  D <- dose[dose.index]
  V<-psi[id[obs.index],1]
  k<-psi[id[obs.index],2]

  V.dose<-psi[id[dose.index],1]
  k.dose<-psi[id[dose.index],2]

  cumsum <- D*exp(-k.dose*t.dose)

  N <- length(obs.index)
  Nsubj <- length(unique(id))
  ypred <- rep(0,N)

  browser()
  i = 1 #obs
  j = 2 #obs per indiv

  for (i in 1:N){
    indiv <- id[dose.index][i]
    for (j in 1:length(obs.index[id==indiv])){
      ypred[i] <- sum(cumsum[id==indiv][0:(obs.index[id==indiv][j])])
    }  
  }
  
  

  
  return(ypred)
}

# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="tobra",type="structural",
  ,psi0=matrix(c(1,7),ncol=2,byrow=TRUE, dimnames=list(NULL, c("V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))


K1 = 100
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2


replicate = 3
seed0 = 395246
seed1 = 3952


final_rwm <- 0
final_mix <- 0
final_mala <- 0


m=1
for (m in 1:replicate){
  print(m)
  l = list(c(50,2),c(80,3),c(60,4))
  # l = list(c(1,5,2,0,0,0),c(3,12,5,0,0,0),c(6,3,7,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model<-saemixModel(model=model1cpt,description="tobra",type="structural"
  ,psi0=matrix(l[[m]],ncol=2,byrow=TRUE, dimnames=list(NULL, c("V","k"))),
  transform.par=c(1,1),omega.init=matrix(c(1/m,0,0,1/m),ncol=2,byrow=TRUE))

  options<-list(seed=seed0/m,map=F,fim=F,ll.is=T,nb.chains = 1,
   nbiter.mcmc = c(2,2,2,0,0,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0)
  tobra_ref<-data.frame(saemix(saemix.model,saemix.data,options)$par)
  tobra_ref <- cbind(iterations, tobra_ref)
  tobra_ref[,4:6] <- sqrt(tobra_ref[,4:6])
  tobra_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,tobra_ref[-1,])

  options.new<-list(seed=m*seed1,map=F,fim=F,ll.is=T,nb.chains = 1,
   nbiter.mcmc = c(2,2,2,2,0,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),
   map.range=c(1:4),nbiter.burn =0)
  tobra_new_ref<-data.frame(saemix(saemix.model,saemix.data,options.new)$par)
  tobra_mix <- cbind(iterations, tobra_new_ref)
  tobra_mix[,4:6] <- sqrt(tobra_mix[,4:6])
  tobra_mix['individual'] <- m
  final_mix <- rbind(final_mix,tobra_mix[-1,])

  #  options.mala<-list(seed=seed0/m,map=F,fim=F,ll.is=T,nb.chains = 1,
  #   nbiter.mcmc = c(2,2,2,0,2,0),
  #   nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(1),nbiter.burn =0,sigma.val=0.002,gamma.val=0.1)
  # tobra_mala_ref<-data.frame(saemix(saemix.model,saemix.data,options.mala)$par)
  # tobra_mala <- cbind(iterations, tobra_mala_ref)
  # tobra_mala[,4:6] <- sqrt(tobra_mala[,4:6])
  # tobra_mala['individual'] <- m
  # final_mala <- rbind(final_mala,tobra_mala[-1,])

  #  options.nuts<-list(seed=seed0/m,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,0,6),
  #   nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(1:15),nbiter.burn =0,sigma.val=0.002,gamma.val=0.1, modelstan = modelstan)
  # tobra_nuts_ref<-data.frame(saemix(saemix.model,saemix.data,options.nuts)$par)
  # tobra_nuts <- cbind(iterations, tobra_nuts_ref)
  # tobra_nuts[,5:7] <- sqrt(tobra_nuts[,5:7])
  # tobra_nuts['individual'] <- m
  # final_nuts <- rbind(final_nuts,tobra_nuts[-1,])

}


convpop <- graphConvMC_diffpk1(final_rwm[,c(1,3,7)],final_mix[,c(1,3,7)])
convvar <- graphConvMC_diffpk1(final_rwm[,c(1,5,7)],final_mix[,c(1,5,7)])

convpop <- graphConvMC_diffpk1(final_rwm[,c(1,2,7)],final_mix[,c(1,2,7)])
convvar <- graphConvMC_diffpk1(final_rwm[,c(1,4,7)],final_mix[,c(1,4,7)])







convpop <- graphConvMC_diffpk2_3df(final_rwm[,c(1,3,9)],final_mix[,c(1,3,9)],final_mala[,c(1,3,9)])
convvar <- graphConvMC_diffpk1_3df(final_rwm[,c(1,6,9)],final_mix[,c(1,6,9)],final_mala[,c(1,6,9)])

convpop <- graphConvMC_diffpk2(final_rwm[,c(1,3,9)],final_mix[,c(1,3,9)])
convvar <- graphConvMC_diffpk1(final_rwm[,c(1,6,9)],final_mix[,c(1,6,9)])

save <- grid.arrange(convpop,convvar, ncol=2)
# ggsave(save,file="pics_square/convpkseednew.pdf", width = 900, height = 225, units = "mm")

# ggsave(save,file="newpics/convpk.pdf", width = 900, height = 450, units = "mm")

ggsave(save,file="/Users/karimimohammedbelhal/Desktop/convpk_mala.pdf", width = 900, height = 250, units = "mm")

convpop2 <- graphConvMC_diffpk2(final_rwm[,c(1,2,9)],final_mix[,c(1,2,9)])
convvar2 <- graphConvMC_diffpk1(final_rwm[,c(1,5,9)],final_mix[,c(1,5,9)])


convpop3 <- graphConvMC_diffpk2(final_rwm[,c(1,4,9)],final_mix[,c(1,4,9)])
convvar3 <- graphConvMC_diffpk1(final_rwm[,c(1,7,9)],final_mix[,c(1,7,9)])


ka_true <- 8
V_true <- 10
k_true <- 0.1
o_ka <- 0.5
o_V <- 0.2
o_k <- 0.2
a_true <- 1

var_rwm <- 0
error_rwm <- 0


var_mix <- 0
error_mix <- 0

replicate = 4

final_rwm <- 0
final_mix <- 0

replicate <- 40
for (m in 1:replicate){
  
  
myModel <- inlineModel("


[INDIVIDUAL]
input = {ka_pop, V_pop, k_pop, omega_ka, omega_V, omega_k}
DEFINITION:
ka = {distribution=lognormal, reference=ka_pop, sd=omega_ka}
V  = {distribution=lognormal, reference=V_pop,  sd=omega_V }
k = {distribution=lognormal, reference=k_pop, sd=omega_k}


[LONGITUDINAL]
input = {ka, V, k,a}
EQUATION:
C = pkmodel(ka,V,k)
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=50

pop.param   <- c(
  ka_pop  = ka_true,    omega_ka  = o_ka,
  V_pop   = V_true,   omega_V   = o_V,
  k_pop  = k_true,    omega_k  = o_k, a =a_true)
  

res <- simulx(model     = myModel,
              parameter = pop.param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(0,10,by=1)))
  

  tobra.saemix <- res$y
  tobra.saemix$amount <- 100
  
  saemix.model<-saemixModel(model=model1cpt,description="tobra",type="structural"
  ,psi0=matrix(c(1,10,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE))


saemix.data<-saemixData(name.data=tobra.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")


  options<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0)
  tobra_ref<-data.frame(saemix(saemix.model,saemix.data,options)$par)
  tobra_ref <- cbind(iterations, tobra_ref)
  tobra_ref[,5:8] <- sqrt(tobra_ref[,5:8])
  ML <- tobra_ref[,2:8]
  ML[1:(end+1),]<- tobra_ref[end+1,2:8]
  error_rwm <- error_rwm + (tobra_ref[,2:8]-ML)^2
  tobra_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,tobra_ref)
  

  options.new<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),map.range=c(1:6),nbiter.burn =0)
  tobra_mix<-data.frame(saemix(saemix.model,saemix.data,options.new)$par)
  tobra_mix <- cbind(iterations, tobra_mix)
  tobra_mix[,5:8] <- sqrt(tobra_mix[,5:8])
  ML <- tobra_mix[,2:8]
  ML[1:(end+1),]<- tobra_mix[end+1,2:8]
  error_mix <- error_mix + (tobra_mix[,2:8]-ML)^2
  tobra_mix['individual'] <- m
  final_mix <- rbind(final_mix,tobra_mix)
}

error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix

error_rwm <- cbind(iterations, error_rwm)
error_mix <- cbind(iterations, error_mix)
error_rwm[2,] = error_mix[2,]


err_mix<- tobra_ref
err_rwm<- tobra_ref
err_rwm[,2:8] <- error_rwm[,2:8]
err_mix[,2:8] <- error_mix[,2:8]


a <- err_rwm[-1,]
b <- err_mix[-1,]


graphConvMC_se1 <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
      xlab("iteration") +scale_x_log10(breaks= c(10,100,200))+ ylab(expression(paste(E(V[pop]))))  
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


graphConvMC_se2 <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
      xlab("iteration") +scale_x_log10(breaks= c(10,100,200))+ ylab(expression(paste(E(omega[V]))))  
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


c <- graphConvMC_se1(a[,c(1,3,9)],b[,c(1,3,9)])
d <- graphConvMC_se2(a[,c(1,6,9)],b[,c(1,6,9)])


save <- grid.arrange(c,d, ncol=2)
# ggsave(save,file="pics_square/sepk.pdf", width = 900, height = 225, units = "mm")

ggsave(save,file="newpics/se_pk.pdf", width = 900, height = 450, units = "mm")
ggsave(save,file="/Users/karimimohammedbelhal/Desktop/se_pk.pdf", width = 900, height = 250, units = "mm")
# c <- graphConvMC_se1(a[,c(1,2,9)],b[,c(1,2,9)])
# d <- graphConvMC_se2(a[,c(1,5,9)],b[,c(1,5,9)])

# grid.arrange(c,d, ncol=2)

# c <- graphConvMC_se1(a[,c(1,4,9)],b[,c(1,4,9)])
# d <- graphConvMC_se2(a[,c(1,7,9)],b[,c(1,7,9)])

# grid.arrange(c,d, ncol=2)


# graphConvMC_diff4(final_rwm[,c(1,3,9)],final_mix[,c(1,3,9)])