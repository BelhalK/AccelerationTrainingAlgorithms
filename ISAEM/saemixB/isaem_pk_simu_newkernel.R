load("minibatch_pk_simu.RData")
# save.image("minibatch_pk_simu.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/R")
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
  
  source('main_incremental.R')
  source('main_estep_incremental.R')


  source('/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/R/mixtureFunctions.R')
  source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/plots.R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB")


library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)


# warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/csda_new/data/warfarin_data.txt", header=T)



model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  k <- Cl/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*time)-exp(-ka*time))
  return(ypred)
}

model <- inlineModel("


[INDIVIDUAL]
input = {ka_pop, V_pop, Cl_pop, omega_ka, omega_V, omega_Cl}
DEFINITION:
ka = {distribution=lognormal, reference=ka_pop, sd=omega_ka}
V  = {distribution=lognormal, reference=V_pop,  sd=omega_V }
Cl = {distribution=lognormal, reference=Cl_pop, sd=omega_Cl}


[LONGITUDINAL]
input = {ka, V, Cl,a}
EQUATION:
C = pkmodel(ka,V,Cl)
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=200

param   <- c(
  ka_pop  = 1,    omega_ka  = 0.3,
  V_pop   = 10,   omega_V   = 0.2,
  Cl_pop  = 1,    omega_Cl  = 0.3, a =1)
  
res <- simulx(model     = model,
              parameter = param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(1,10,by=1)))

 warfarin.saemix <- res$y
 warfarin.saemix$amount <- 100
 saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")

# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(3,3,0.1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),fixed.estim=c(1,1,1),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))



K1 = 300
K2 = 100
iterations = 0:(K1+K2-1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50

seed0=3456


options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='seq')
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])


# options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
#   nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
#   nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
# theo_mix75<-saemix_incremental(saemix.model,saemix.data,options.incremental75)
# theo_mix75 <- data.frame(theo_mix75$param)
# theo_mix75 <- cbind(iterations, theo_mix75[-1,])


# options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), 
#                           nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,
#                           nbiter.burn =0, nb.replacement=50,sampling='randompass')
# theo_mix50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
# theo_mix50 <- data.frame(theo_mix50$param)
# theo_mix50 <- cbind(iterations, theo_mix50[-1,])


# options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
#   nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
#   nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
# theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
# theo_mix25 <- data.frame(theo_mix25$param)
# theo_mix25 <- cbind(iterations, theo_mix25[-1,])

# theo_ref_scaled <- theo_ref
# theo_mix25_scaled <- theo_mix25
# theo_mix50_scaled <- theo_mix50
# theo_mix75_scaled <- theo_mix75
# theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
# theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25
# theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
# theo_mix75_scaled$iterations = theo_mix75_scaled$iterations*0.75

# graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled)

# graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled)
# graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix50_scaled,theo_mix75_scaled)

###NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL####

options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,2),
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, 
 map.range=c(1:4), nb.replacement=100,sampling='seq')
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])



options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,2), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:4),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,2), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(1:4),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])


theo_ref_scaled <- theo_ref
theo_mix50_scaled <- theo_mix50
theo_mix25_scaled <- theo_mix25
theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25
graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled)

graphConvMC_threekernels(theo_ref,theo_mix50,theo_mix25)
graphConvMC_threekernels(theo_ref_scaled,theo_mix25_scaled,theo_mix25_scaled)
###NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL#######NEWKERNEL####


options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
theo_mix75<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental75))
theo_mix75 <- cbind(iterations, theo_mix75[-1,])


options.incremental85<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=85,sampling='randompass')
theo_mix85<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental85))
theo_mix85 <- cbind(iterations, theo_mix85[-1,])


theo_ref_scaled <- theo_ref
theo_mix50_scaled <- theo_mix50
theo_mix25_scaled <- theo_mix25
theo_mix75_scaled <- theo_mix75
theo_mix85_scaled <- theo_mix85


theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25
theo_mix75_scaled$iterations = theo_mix75_scaled$iterations*0.75
theo_mix85_scaled$iterations = theo_mix85_scaled$iterations*0.85

# graphConvMC_threekernels(theo_ref_scaled,theo_mix50_scaled,theo_mix50_scaled)
# graphConvMC_threekernels(theo_ref_scaled,theo_mix50_scaled,theo_mix25_scaled)
graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix75_scaled,theo_mix85_scaled)



final_ref <- 0
final_mix <- 0

ka_true <- 1
V_true <- 10
Cl_true <- 1
a_true <- 1

o_ka <- 0.3 #o^2=0.09
o_V <- 0.2  #o^2=0.04
o_Cl <- 0.3  #o^2=0.09

error_rwm <- 0
error_mix50 <- 0
error_mix25 <- 0
error_mix75 <- 0

K1 = 200
K2 = 100
iterations = 0:(K1+K2-1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50


map_range = 1:4
true_param <- data.frame("ka" = ka_true, "V" = V_true, "Cl" = Cl_true, "omega2.ka"=o_ka^2 ,"omega2.V"= o_V^2,"omega2.Cl"= o_Cl^2, "a" = a_true)
seed0 = 39546
replicate = 20
for (m in 1:replicate){

model <- inlineModel("


[INDIVIDUAL]
input = {ka_pop, V_pop, Cl_pop, omega_ka, omega_V, omega_Cl}
DEFINITION:
ka = {distribution=lognormal, reference=ka_pop, sd=omega_ka}
V  = {distribution=lognormal, reference=V_pop,  sd=omega_V }
Cl = {distribution=lognormal, reference=Cl_pop, sd=omega_Cl}


[LONGITUDINAL]
input = {ka, V, Cl,a}
EQUATION:
C = pkmodel(ka,V,Cl)
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=200

param   <- c(
  ka_pop  = ka_true,    omega_ka  = o_ka,
  V_pop   = V_true,   omega_V   = o_V,
  Cl_pop  = Cl_true,    omega_Cl  = o_Cl, a =a_true)
  
res <- simulx(model     = model,
              parameter = param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(1,10,by=1)))

 warfarin.saemix <- res$y
 warfarin.saemix$amount <- 100
 saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")

# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(3,3,0.1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),fixed.estim=c(1,1,1),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


options<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,2),
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, 
 map.range=c(map_range), nb.replacement=100,sampling='seq', gamma= 0)
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])
ML <- theo_ref[,2:8]
ML[0:end,1:7]<- theo_ref[end,2:8]
# ML[0:end,1:7]<- true_param
error_rwm <- error_rwm + (theo_ref[,2:8]-ML)^2
theo_ref['individual'] <- m


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,2), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(map_range),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])

ML <- theo_ref[,2:8]
ML[0:end,1:7]<- theo_mix50[end,2:8]
# ML[0:end,1:7]<- true_param
error_mix50 <- error_mix50 + (theo_mix50[,2:8]-ML)^2
theo_mix50['individual'] <- m


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,2), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(map_range),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])


ML <- theo_ref[,2:8]
ML[0:end,1:7]<- theo_mix25[end,2:8]
# ML[0:end,1:7]<- true_param
error_mix25 <- error_mix25 + (theo_mix25[,2:8]-ML)^2
theo_mix25['individual'] <- m


options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,2), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(map_range),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
theo_mix75<-saemix_incremental(saemix.model,saemix.data,options.incremental75)
theo_mix75 <- data.frame(theo_mix75$param)
theo_mix75 <- cbind(iterations, theo_mix75[-1,])


ML <- theo_ref[,2:8]
ML[0:end,1:7]<- theo_mix75[end,2:8]
# ML[0:end,1:7]<- true_param
error_mix75 <- error_mix75 + (theo_mix75[,2:8]-ML)^2
theo_mix75['individual'] <- m


}


error_rwm <- 1/replicate*error_rwm
error_mix25 <- 1/replicate*error_mix25
error_mix50 <- 1/replicate*error_mix50
error_mix75 <- 1/replicate*error_mix75

err_rwm<- theo_ref
err_mix25<- theo_ref
err_mix50<- theo_ref
err_mix75<- theo_ref

err_rwm[,2:8] <- error_rwm
err_mix25[,2:8] <- error_mix25
err_mix50[,2:8] <- error_mix50
err_mix75[,2:8] <- error_mix75


err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = err_rwm_scaled$iterations*1


err_mix25_scaled <- err_mix25
err_mix25_scaled$iterations = err_mix25_scaled$iterations*0.25


err_mix50_scaled <- err_mix50
err_mix50_scaled$iterations = err_mix50_scaled$iterations*0.5

err_mix75_scaled <- err_mix75
err_mix75_scaled$iterations = err_mix75_scaled$iterations*0.75

err_rwm_scaled$algo <- 'SAEM'
err_mix25_scaled$algo <- 'ISAEM25'
err_mix50_scaled$algo <- 'ISAEM50'
err_mix75_scaled$algo <- 'ISAEM75'

graphConvMC_se1 <- function(df,df2, df3,df4,title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  df3$individual <- as.factor(df3$individual)
  df4$individual <- as.factor(df4$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue",linetype = 1,size=0.8)+
    geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red",linetype = 1,size=0.8)+geom_line(aes_string(df4[,1],df4[,j],by=df4[,ncol(df4)]),colour="yellow",linetype = 1,size=0.8)+
      xlab("iteration") +scale_x_log10()+ ylab(expression(paste(E(V[pop]))))  
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


c <- graphConvMC_se1(err_rwm_scaled[,c(1,2,9)],err_mix25_scaled[,c(1,2,9)],err_mix50_scaled[,c(1,2,9)],err_mix75_scaled[,c(1,2,9)])
c <- graphConvMC_se1(err_rwm_scaled[,c(1,3,9)],err_mix25_scaled[,c(1,3,9)],err_mix50_scaled[,c(1,3,9)],err_mix75_scaled[,c(1,3,9)])
d <- graphConvMC_se1(err_rwm_scaled[,c(1,6,9)],err_mix25_scaled[,c(1,6,9)],err_mix50_scaled[,c(1,6,9)],err_mix75_scaled[,c(1,6,9)])


save <- grid.arrange(c,d, ncol=2)