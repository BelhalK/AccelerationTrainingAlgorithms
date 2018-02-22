
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

warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")



model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
  Tlag<-psi[id,1]
  ka<-psi[id,2]
  V<-psi[id,3]
  Cl<-psi[id,4]
  k<-Cl/V
  dt <- pmax(time-Tlag, 0)
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*dt)-exp(-ka*dt))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.2,1,7,1),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","Cl"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE))




K1 = 600
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50

final_rwm <- 0
final_ref <- 0
var_rwm <- 0
final_mix <- 0
final_mix25 <- 0
var_mix <- 0
var_mix25 <- 0

error_rwm <- 0
error_mixseq <- 0
error_mix25seq <- 0

error_mixiter <- 0
error_mix25iter <- 0

error_mixpass <- 0
error_mix25pass <- 0

Tlag_true=0.78
ka_true <- 1
V_true <- 8
Cl_true <- 0.1

o_Tlag <- 0.57 #o^2=0.32
o_ka <- 0.5 #o^2=0.25
o_V <- 0.2  #o^2=0.04
o_Cl <- 0.3  #o^2=0.09
a_true = 0.266
beta_Cl_lw70_true = 0.60411
beta_V_lw70_true = 0.8818

true_param <- data.frame("Tlag" = Tlag_true,"ka" = ka_true, "V" = V_true, "Cl" = Cl_true, "omega2.Tlag"=o_Tlag^2, "omega2.ka"=o_ka^2 ,"omega2.V"= o_V^2,"omega2.Cl"= o_Cl^2, "a" = a_true)
seed0 = 39546
replicate = 20
for (m in 1:replicate){
  
  myModel <- inlineModel("

[COVARIATE]
input = wt

EQUATION:
lw70 = log(wt/70)

[INDIVIDUAL]
input = {Tlag_pop, omega_Tlag, ka_pop, omega_ka, V_pop, beta_V_lw70, lw70, omega_V, Cl_pop, beta_Cl_lw70, omega_Cl}

DEFINITION:
Tlag = {distribution=lognormal, typical=Tlag_pop, sd=omega_Tlag}
ka = {distribution=lognormal, typical=ka_pop, sd=omega_ka}
V = {distribution=lognormal, typical=V_pop, covariate=lw70, coefficient=beta_V_lw70, sd=omega_V}
Cl = {distribution=lognormal, typical=Cl_pop, covariate=lw70, coefficient=beta_Cl_lw70, sd=omega_Cl}

[LONGITUDINAL]
input =  {Tlag, ka, V, Cl,a}

EQUATION:
Cc = pkmodel(Tlag, ka, V, Cl)

OUTPUT:
output = {Cc}

DEFINITION:
y1 = {distribution=normal, prediction=Cc, sd=a}
")


populationParameter   <- c(Tlag_pop= Tlag_true, omega_Tlag= o_Tlag,
  ka_pop  = ka_true,    omega_ka  = o_ka,
  V_pop   = V_true,   omega_V   = o_V,
  Cl_pop  = Cl_true,    omega_Cl  = o_Cl, a =a_true, beta_V_lw70 = beta_V_lw70_true, beta_Cl_lw70 = beta_Cl_lw70_true)


trt <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/treatment.txt", header = TRUE) 
originalId<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/originalId.txt', header=TRUE) 
individualCovariate<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/individualCovariate.txt', header = TRUE) 
time<-read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/output1.txt",header=TRUE)


list.param <- list(populationParameter,individualCovariate)
name<-"y1"
out1<-list(name=name,time=time) 

# call the simulator 
res <- simulx(model=myModel,treatment=trt,parameter=list.param,output=out1)
individualCovariate$wt <- log(individualCovariate$wt/70)
warfarin.saemix <- res$y1
treat <- res$treatment[,c(1,3)]
covandtreat <- merge(individualCovariate ,treat,by="id")
warfarin.saemix <- merge(covandtreat ,warfarin.saemix,by="id")


saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.2,3,10,2),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","Cl"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE),covariate.model=matrix(c(0,0,1,1),ncol=4,byrow=TRUE),error.model="constant")

saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time", name.covariates=c("wt"),units=list(x="kg",
  covariates=c("kg/ha")))



options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='')
theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

ML <- theo_ref[,2:10]
# ML[1:(end+1),]<- theo_ref[end+1,2:10]
ML[1:(end+1),1:9]<- true_param
error_rwm <- error_rwm + (theo_ref[,2:10]-ML)^2
theo_ref['individual'] <- m
final_ref <- rbind(final_ref,theo_ref)

#SEQ
options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='seq')
theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
theo_mix <- cbind(iterations, theo_mix)
ML <- theo_mix[,2:10]
# ML[1:(end+1),]<- theo_mix[end+1,2:10]
ML[1:(end+1),1:9]<- true_param
error_mixseq <- error_mixseq + (theo_mix[,2:10]-ML)^2
theo_mix['individual'] <- m
final_mix <- rbind(final_mix,theo_mix)

options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='seq')
theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_mix25)
ML <- theo_mix25[,2:10]
# ML[1:(end+1),]<- theo_mix25[end+1,2:10]
ML[1:(end+1),1:9]<- true_param
error_mix25seq <- error_mix25seq + (theo_mix25[,2:10]-ML)^2
theo_mix25['individual'] <- m
final_mix25 <- rbind(final_mix25,theo_mix25)

# #RANDOMPASS
# options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randompass')
# theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
# theo_mix <- cbind(iterations, theo_mix)
# ML <- theo_mix[,2:10]
# # ML[1:(end+1),]<- theo_mix[end+1,2:10]
# ML[1:(end+1),1:9]<- true_param
# error_mixpass <- error_mixpass + (theo_mix[,2:10]-ML)^2
# theo_mix['individual'] <- m
# final_mix <- rbind(final_mix,theo_mix)

# options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randompass')
# theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))  
# theo_mix25 <- cbind(iterations, theo_mix25)
# ML <- theo_mix25[,2:10]
# # ML[1:(end+1),]<- theo_mix25[end+1,2:10]
# ML[1:(end+1),1:9]<- true_param
# error_mix25pass <- error_mix25pass + (theo_mix25[,2:10]-ML)^2
# theo_mix25['individual'] <- m
# final_mix25 <- rbind(final_mix25,theo_mix25)

# #RANDOMITER
# options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randomiter')
# theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
# theo_mix <- cbind(iterations, theo_mix)
# ML <- theo_mix[,2:10]
# # ML[1:(end+1),]<- theo_mix[end+1,2:10]
# ML[1:(end+1),1:9]<- true_param
# error_mixiter <- error_mixiter + (theo_mix[,2:10]-ML)^2
# theo_mix['individual'] <- m
# final_mix <- rbind(final_mix,theo_mix)

# options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randomiter')
# theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
# theo_mix25 <- cbind(iterations, theo_mix25)
# ML <- theo_mix25[,2:10]
# # ML[1:(end+1),]<- theo_mix25[end+1,2:10]
# ML[1:(end+1),1:9]<- true_param
# error_mix25iter <- error_mix25iter + (theo_mix25[,2:10]-ML)^2
# theo_mix25['individual'] <- m
# final_mix25 <- rbind(final_mix25,theo_mix25)
 
}

graphConvMC_diff(final_ref,final_ref,final_ref)
# graphConvMC_diff(final_ref,final_mix,final_mix25)

K_inc <- 12
error_rwm <- 1/replicate*error_rwm
err_rwm<- theo_ref[-1,]
err_rwm[,2:10] <- error_rwm[-1,]
err_rwm_scaled <- err_rwm
err_rwm_scaled[1:(K_inc/4),]$iterations = seq(1, K_inc, by=4)
err_rwm_scaled[((K_inc/4)+1):end,]$iterations = err_rwm_scaled[((K_inc/4)+1):end,]$iterations + err_rwm_scaled[(K_inc/4),]$iterations 
err_rwm_scaled$algo <- 'SAEM'
err_rwm_scaled$method <- 'seq'



error_mixseq <- 1/replicate*error_mixseq
err_mixseq<- theo_ref[-1,]
err_mixseq[,2:10] <- error_mixseq[-1,]
err_mixseq_scaled <- err_mixseq
err_mixseq_scaled[1:(K_inc/2),]$iterations = seq(1, K_inc, by=2)
err_mixseq_scaled[((K_inc/2)+1):end,]$iterations = err_mixseq_scaled[((K_inc/2)+1):end,]$iterations + err_mixseq_scaled[(K_inc/2),]$iterations 
err_mixseq_scaled$algo <- 'ISAEM50'
err_mixseq_scaled$method <- 'seq'


error_mixpass <- 1/replicate*error_mixpass
err_mixpass<- theo_ref[-1,]
err_mixpass[,2:10] <- error_mixpass[-1,]
err_mixpass_scaled <- err_mixpass
err_mixpass_scaled[1:(K_inc/2),]$iterations = seq(1, K_inc, by=2)
err_mixpass_scaled[((K_inc/2)+1):end,]$iterations = err_mixpass_scaled[((K_inc/2)+1):end,]$iterations + err_mixpass_scaled[(K_inc/2),]$iterations 
err_mixpass_scaled$algo <- 'ISAEM50'
err_mixpass_scaled$method <- 'pass'


error_mixiter <- 1/replicate*error_mixiter
err_mixiter<- theo_ref[-1,]
err_mixiter[,2:10] <- error_mixiter[-1,]
err_mixiter_scaled <- err_mixiter
err_mixiter_scaled[1:(K_inc/2),]$iterations = seq(1, K_inc, by=2)
err_mixiter_scaled[((K_inc/2)+1):end,]$iterations = err_mixiter_scaled[((K_inc/2)+1):end,]$iterations + err_mixiter_scaled[(K_inc/2),]$iterations 
err_mixiter_scaled$algo <- 'ISAEM50'
err_mixiter_scaled$method <- 'iter'



error_mix25seq <- 1/replicate*error_mix25seq
err_mix25seq<- theo_ref[-1,]
err_mix25seq[,2:10] <- error_mix25seq[-1,]
err_mix25seq$iterations = 1:((K1+K2))
err_mix25seq$algo <- 'ISAEM25'
err_mix25seq$method <- 'seq'



error_mix25pass <- 1/replicate*error_mix25pass
err_mix25pass<- theo_ref[-1,]
err_mix25pass[,2:10] <- error_mix25pass[-1,]
err_mix25pass$iterations = 1:((K1+K2))
err_mix25pass$algo <- 'ISAEM25'
err_mix25pass$method <- 'pass'



error_mix25iter <- 1/replicate*error_mix25iter
err_mix25iter<- theo_ref[-1,]
err_mix25iter[,2:10] <- error_mix25iter[-1,]
err_mix25iter$iterations = 1:((K1+K2))
err_mix25iter$algo <- 'ISAEM25'
err_mix25iter$method <- 'iter'



for (i in 2:10){
# i = 6
comparison <- 0
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,12,13)],err_mixseq_scaled [0:end,c(1,i,12,13)],err_mix25seq[0:end,c(1,i,12,13)],
                                              err_mixpass_scaled [0:end,c(1,i,12,13)],err_mix25pass[0:end,c(1,i,12,13)],
                                              err_mixiter_scaled [0:end,c(1,i,12,13)],err_mix25iter[0:end,c(1,i,12,13)])

var <- melt(comparison, id.var = c('iterations','algo','method'), na.rm = TRUE)


prec <- seplot(var, title="ALGO - EM (same complexity)",legend=TRUE)
setwd("/Users/karimimohammedbelhal/Desktop/")
ggsave(paste("precwarfa_", i, ".png", sep=""),prec)
}

error_rwm <- 1/replicate*error_rwm
err_rwm<- theo_ref[-1,]
err_rwm[,2:10] <- error_rwm[-1,]
err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = seq(1, 4*end, by=4)
err_rwm_scaled$algo <- 'SAEM'
err_rwm_scaled$method <- 'seq'

error_mixseq <- 1/replicate*error_mixseq
err_mixseq<- theo_ref[-1,]
err_mixseq[,2:10] <- error_mixseq[-1,]
err_mixseq_scaled <- err_mixseq
err_mixseq_scaled$iterations = seq(1, 2*end, by=2)
err_mixseq_scaled$algo <- 'ISAEM50'
err_mixseq_scaled$method <- 'seq'

error_mix25seq <- 1/replicate*error_mix25seq
err_mix25seq<- theo_ref[-1,]
err_mix25seq[,2:10] <- error_mix25seq[-1,]
err_mix25seq$iterations = 1:((K1+K2))
err_mix25seq$algo <- 'ISAEM25'
err_mix25seq$method <- 'seq'


for (i in 2:10){
# i = 6
comparison <- 0
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,12,13)],err_mixseq_scaled [0:end,c(1,i,12,13)],err_mix25seq[0:end,c(1,i,12,13)])

var <- melt(comparison, id.var = c('iterations','algo','method'), na.rm = TRUE)


prec <- seplot(var, title="ALGO - EM (same complexity)",legend=TRUE)
# assign(paste("prec", i, sep = ""), prec) 
# setwd("/Users/karimimohammedbelhal/Desktop/")
# ggsave(paste("precwarfa_seq_50sim_100indiv_", i, ".png", sep=""),prec)
}

grid.arrange(prec2,prec3,prec4,prec5,prec6,prec7,prec8,prec9, ncol=4)
