
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
  k <- Cl/V
  dt <- pmax(time-Tlag, 0)
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*dt)-exp(-ka*dt))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.2,1,7,1),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","Cl"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE),error.model="combined")



K1 = 300
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


populationParameter<- read.vector('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/populationParameter2.txt') 
true_param <- populationParameter[c(1,2,3,5,11,12,13,14,19,20)]
true_param <- data.frame("Tlag" = true_param[1],
                          "ka" = true_param[2],
                          "V" = true_param[3], 
                          "Cl" = true_param[4], 
                          "omega2.Tlag"=true_param[5]^2, 
                          "omega2.ka"=true_param[6]^2 ,
                          "omega2.V"= true_param[7]^2,"omega2.Cl"= true_param[8]^2, 
                          "a" = true_param[9],"b"=true_param[10])
seed0 = 39546
replicate = 50
for (m in 1:replicate){
  

      myModel <- inlineModel("[COVARIATE]
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
input = {a1, b1}
 
DESCRIPTION: PK oral + indirect response model

input =  {Tlag, ka, V, Cl}

EQUATION:
Cc = pkmodel(Tlag, ka, V, Cl)

OUTPUT:
output = {Cc}

DEFINITION:
y1 = {distribution=normal, prediction=Cc, errorModel=combined1(a1, b1)}

")
# treatment
trt <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/treatment.txt", header = TRUE) 

# parameters 
originalId<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/originalId.txt', header=TRUE) 
populationParameter<- read.vector('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/populationParameter2.txt') 
individualCovariate<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/individualCovariate.txt', header = TRUE) 
list.param <- list(populationParameter,individualCovariate)
# output 
name<-"y1"
time<-read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/output1.txt",header=TRUE)
out1<-list(name=name,time=time) 



# call the simulator 
res <- simulx(model=myModel,treatment=trt,parameter=list.param,output=out1)

individualCovariate[1:10,]$id <- 1:10
warfarin.saemix <- res$y1
treat <- res$treatment[,c(1,3)]
covandtreat <- merge(individualCovariate ,treat,by="id")
warfarin.saemix <- merge(covandtreat ,warfarin.saemix,by="id")


saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.2,3,13,2),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","Cl"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE),covariate.model=matrix(c(0,0,1,1),ncol=4,byrow=TRUE),error.model="combined")

saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time", name.covariates=c("wt"))


options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='')
theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

ML <- theo_ref[,2:11]
# ML[1:(end+1),]<- theo_ref[end+1,2:10]
ML[1:(end+1),1:10]<- true_param
error_rwm <- error_rwm + (theo_ref[,2:11]-ML)^2
theo_ref['individual'] <- m
final_ref <- rbind(final_ref,theo_ref)

# #SEQ
# options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='seq')
# theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
# theo_mix <- cbind(iterations, theo_mix)
# ML <- theo_mix[,2:11]
# # ML[1:(end+1),]<- theo_mix[end+1,2:11]
# ML[1:(end+1),1:10]<- true_param
# error_mixseq <- error_mixseq + (theo_mix[,2:11]-ML)^2
# theo_mix['individual'] <- m
# final_mix <- rbind(final_mix,theo_mix)

# options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='seq')
# theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
# theo_mix25 <- cbind(iterations, theo_mix25)
# ML <- theo_mix25[,2:11]
# # ML[1:(end+1),]<- theo_mix25[end+1,2:11]
# ML[1:(end+1),1:10]<- true_param
# error_mix25seq <- error_mix25seq + (theo_mix25[,2:11]-ML)^2
# theo_mix25['individual'] <- m
# final_mix25 <- rbind(final_mix25,theo_mix25)

#RANDOMPASS
options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randompass')
theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
theo_mix <- cbind(iterations, theo_mix)
ML <- theo_mix[,2:11]
# ML[1:(end+1),]<- theo_mix[end+1,2:11]
ML[1:(end+1),1:10]<- true_param
error_mixpass <- error_mixpass + (theo_mix[,2:11]-ML)^2
theo_mix['individual'] <- m
final_mix <- rbind(final_mix,theo_mix)

# options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randompass')
# theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))  
# theo_mix25 <- cbind(iterations, theo_mix25)
# ML <- theo_mix25[,2:11]
# # ML[1:(end+1),]<- theo_mix25[end+1,2:11]
# ML[1:(end+1),1:10]<- true_param
# error_mix25pass <- error_mix25pass + (theo_mix25[,2:11]-ML)^2
# theo_mix25['individual'] <- m
# final_mix25 <- rbind(final_mix25,theo_mix25)

# #RANDOMITER
# options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randomiter')
# theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
# theo_mix <- cbind(iterations, theo_mix)
# ML <- theo_mix[,2:11]
# # ML[1:(end+1),]<- theo_mix[end+1,2:11]
# ML[1:(end+1),1:10]<- true_param
# error_mixiter <- error_mixiter + (theo_mix[,2:11]-ML)^2
# theo_mix['individual'] <- m
# final_mix <- rbind(final_mix,theo_mix)

options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randomiter')
theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_mix25)
ML <- theo_mix25[,2:11]
# ML[1:(end+1),]<- theo_mix25[end+1,2:11]
ML[1:(end+1),1:10]<- true_param
error_mix25iter <- error_mix25iter + (theo_mix25[,2:11]-ML)^2
theo_mix25['individual'] <- m
final_mix25 <- rbind(final_mix25,theo_mix25)
 
}

# graphConvMC_diff(final_ref,final_ref,final_ref)
# graphConvMC_diff(final_ref,final_mix,final_mix25)

error_rwm <- 1/replicate*error_rwm
error_mixseq <- 1/replicate*error_mixseq
error_mix25seq <- 1/replicate*error_mix25seq

error_mixpass <- 1/replicate*error_mixpass
error_mix25pass <- 1/replicate*error_mix25pass

error_mixiter <- 1/replicate*error_mixiter
error_mix25iter <- 1/replicate*error_mix25iter


err_rwm<- theo_ref[-1,]
err_mixseq<- theo_ref[-1,]
err_mix25seq<- theo_ref[-1,]

err_mixpass<- theo_ref[-1,]
err_mix25pass<- theo_ref[-1,]

err_mixiter<- theo_ref[-1,]
err_mix25iter<- theo_ref[-1,]

err_rwm[,2:11] <- error_rwm[-1,]

err_mixseq[,2:11] <- error_mixseq[-1,]
err_mix25seq[,2:11] <- error_mix25seq[-1,]

err_mixpass[,2:11] <- error_mixpass[-1,]
err_mix25pass[,2:11] <- error_mix25pass[-1,]

err_mixiter[,2:11] <- error_mixiter[-1,]
err_mix25iter[,2:11] <- error_mix25iter[-1,]


err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = seq(1, 4*end, by=4)

err_mixseq_scaled <- err_mixseq
err_mixseq_scaled$iterations = seq(1, 2*end, by=2)

err_mixpass_scaled <- err_mixpass
err_mixpass_scaled$iterations = seq(1, 2*end, by=2)

err_mixiter_scaled <- err_mixiter
err_mixiter_scaled$iterations = seq(1, 2*end, by=2)

err_mix25seq$iterations = 1:((K1+K2))
err_mix25pass$iterations = 1:((K1+K2))
err_mix25iter$iterations = 1:((K1+K2))




err_rwm_scaled$algo <- 'SAEM'
err_rwm_scaled$method <- 'seq'

err_mixseq_scaled$algo <- 'ISAEM50'
err_mixpass_scaled$algo <- 'ISAEM50'
err_mixiter_scaled$algo <- 'ISAEM50'

err_mix25seq$algo <- 'ISAEM25'
err_mix25pass$algo <- 'ISAEM25'
err_mix25iter$algo <- 'ISAEM25'

err_mixseq_scaled$method <- 'seq'
err_mixpass_scaled$method <- 'pass'
err_mixiter_scaled$method <- 'iter'

err_mix25seq$method <- 'seq'
err_mix25pass$method <- 'pass'
err_mix25iter$method <- 'iter'



for (i in 2:10){
# i = 6
comparison <- 0
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,13,14)],err_mixseq_scaled [0:end,c(1,i,13,14)],err_mix25seq[0:end,c(1,i,13,14)],
                                              err_mixpass_scaled [0:end,c(1,i,13,14)],err_mix25pass[0:end,c(1,i,13,14)],
                                              err_mixiter_scaled [0:end,c(1,i,13,14)],err_mix25iter[0:end,c(1,i,13,14)])

var <- melt(comparison, id.var = c('iterations','algo','method'), na.rm = TRUE)


prec <- seplot(var, title="ALGO - EM (same complexity)",legend=TRUE)
# setwd("/Users/karimimohammedbelhal/Desktop/")
# ggsave(paste("precwarfa_", i, ".png", sep=""),prec)
}


# c <- graphConvMC_se2(err_rwm_scaled[-1,c(1,2,8)],err_mix_scaled[-1,c(1,2,8)],err_mix25[-1,c(1,2,8)])
# d <- graphConvMC_se2(err_rwm_scaled[-1,c(1,3,8)],err_mix_scaled[-1,c(1,3,8)],err_mix25[-1,c(1,3,8)])

# grid.arrange(c,d, ncol=2)

# c <- graphConvMC_se2(err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)])
c <- graphConvMC_sec(err_rwm_scaled[0:end,c(1,2,9)],err_mix_scaled [0:end,c(1,2,9)],err_mix25[0:end,c(1,2,9)])
d <- graphConvMC_sed(err_rwm_scaled[0:end,c(1,5,9)],err_mix_scaled[0:end,c(1,5,9)],err_mix25[0:end,c(1,5,9)])

grid.arrange(c,d, ncol=2)

e <- graphConvMC_sec(err_rwm_scaled[0:end,c(1,3,9)],err_mix_scaled[0:end,c(1,3,9)],err_mix25[0:end,c(1,3,9)])
f <- graphConvMC_sed(err_rwm_scaled[0:end,c(1,6,9)],err_mix_scaled[0:end,c(1,6,9)],err_mix25[0:end,c(1,6,9)])

grid.arrange(e,f, ncol=2)

g <- graphConvMC_sec(err_rwm_scaled[0:end,c(1,4,9)],err_mix_scaled[0:end,c(1,4,9)],err_mix25[0:end,c(1,4,9)])
h <- graphConvMC_sed(err_rwm_scaled[0:end,c(1,7,9)],err_mix_scaled[0:end,c(1,7,9)],err_mix25[0:end,c(1,7,9)])

grid.arrange(g,h, ncol=2)



error_rwm <- 1/replicate*error_rwm
error_mixpass <- 1/replicate*error_mixpass

err_rwm<- theo_ref[-1,]
err_mixpass<- theo_ref[-1,]

err_rwm[,2:11] <- error_rwm[-1,]
err_mixpass[,2:11] <- error_mixpass[-1,]


err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = seq(1, 4*end, by=4)

err_mixpass_scaled <- err_mixpass
err_mixpass_scaled$iterations = seq(1, 2*end, by=2)




err_rwm_scaled$algo <- 'SAEM'
err_rwm_scaled$method <- 'seq'

err_mixpass_scaled$algo <- 'ISAEM50'
err_mixpass_scaled$method <- 'pass'


error_mix25iter <- 1/replicate*error_mix25iter
err_mix25iter<- theo_ref[-1,]
err_mix25iter[,2:11] <- error_mix25iter[-1,]
err_mix25iter$iterations = 1:((K1+K2))

err_mix25iter$algo <- 'ISAEM25'
err_mix25iter$method <- 'iter'


for (i in 2:10){
# i = 6
comparison <- 0
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,13,14)],err_mixpass_scaled [0:end,c(1,i,13,14)],err_mix25iter[0:end,c(1,i,13,14)])

var <- melt(comparison, id.var = c('iterations','algo','method'), na.rm = TRUE)


prec <- seplot(var, title="ALGO - EM (same complexity)",legend=TRUE)
# setwd("/Users/karimimohammedbelhal/Desktop/")
# ggsave(paste("precwarfa_", i, ".png", sep=""),prec)
}
