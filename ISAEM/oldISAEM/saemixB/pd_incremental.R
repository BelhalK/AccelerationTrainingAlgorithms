
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
  source('mixtureFunctions.R')

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/")
source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/plots.R")


library("rJava")
library("rCMA")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

PD1.saemix<-read.table( "data/PD1.saemix.tab",header=T,na=".")
PD1.saemix <- subset(PD1.saemix, dose!="90")


PD2.saemix<-read.table( "data/PD2.saemix.tab",header=T,na=".")
saemix.data1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))
saemix.data2<-saemixData(name.data=PD2.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))


modelemax<-function(psi,id,xidep) {
# input:
# psi : matrix of parameters (3 columns, E0, Emax, EC50)
# id : vector of indices
# xidep : dependent variables (same nb of rows as length of id)
# returns:
# a vector of predictions of length equal to length of id
dose<-xidep[,1]
e0<-psi[id,1]
emax<-psi[id,2]
e50<-psi[id,3]
f<-e0+emax*dose/(e50+dose)
return(f)
}


saemix.model<-saemixModel(model=modelemax,description="Emax model",type="structural",
psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")



K1 = 500
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize50<-50
batchsize25<-25

#Weibull
options_pd<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
pd<-data.frame(saemix(saemix.model,saemix.data2,options_pd))
pd<-cbind(iterations,pd)

options_pdincr<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
pdincr<-data.frame(saemix_incremental(saemix.model,saemix.data2,options_pdincr))
pdincr<-cbind(iterations,pdincr)

graphConvMC2_saem(pd,pdincr, title="new kernel")


pd$algo <- 'rwm'
pdincr$algo <- 'ISAEM'

pd_scaled <- pd[rep(seq_len(nrow(pd)), each=100/batchsize50),]
pd_scaled$iterations = 1:(2*(K1+K2+1))


comparison <- 0
# comparison <- rbind(theo_ref,theo_incremental)
comparison <- rbind(pd_scaled[iterations,],pdincr)

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)

options_pdnew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5))
pdnew<-data.frame(saemix(saemix.model,saemix.data2,options_pdnew))



replicate = 3
seed0 = 395246

#RWM
final_rwm <- 0
final_50 <- 0
final_25 <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(20,300,20),c(10,270,10),c(25,320,25),c(20,300,20))


saemix.model<-saemixModel(model=modelemax,description="Emax model",type="structural",
psi0=matrix(l[[m]],ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")

  options_pd<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
  pd<-data.frame(saemix(saemix.model,saemix.data2,options_pd))
  pd<-cbind(iterations,pd)
  pd['individual'] <- m
  pd$algo <- 'rwm'
  pd_scaled <- pd[-1,]
  pd_scaled$iterations = seq(1, 4*end, by=4)
  pd_scaled <- pd_scaled[rep(seq_len(nrow(pd_scaled)), each=4),]
  final_rwm <- rbind(final_rwm,pd_scaled[0:end,])

  options_pdincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
  pdincr25<-data.frame(saemix_incremental(saemix.model,saemix.data2,options_pdincr25))
  pdincr25<-cbind(iterations,pdincr25)
  pdincr25['individual'] <- m
  pdincr25$algo <- 'ISAEM25'
  pdincr25 <- pdincr25[-1,]
  pdincr25$iterations = seq(1, end, by=1)
  final_25 <- rbind(final_25,pdincr25[0:end,])

  options_pdincr50<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
  pdincr50<-data.frame(saemix_incremental(saemix.model,saemix.data2,options_pdincr50))
  pdincr50<-cbind(iterations,pdincr50)
  pdincr50['individual'] <- m
  pdincr50$algo <- 'ISAEM50'
  pd_scaled50 <- pdincr50[-1,]
  pd_scaled50$iterations = seq(1, 2*end, by=2)
  pd_scaled50 <- pd_scaled50[rep(seq_len(nrow(pd_scaled50)), each=2),]
  final_50 <- rbind(final_50,pd_scaled50[0:end,])
}


a <- graphConvMC_diffz(final_rwm[,c(1,3,9)],final_50[,c(1,3,9)],final_25[,c(1,3,9)])
b <- graphConvMC_diffw(final_rwm[,c(1,6,9)],final_50[,c(1,6,9)],final_25[,c(1,6,9)])

grid.arrange(a,b, ncol=2)



K1 = 300
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize50<-50
batchsize25<-25


final_rwm <- 0
final_ref <- 0
var_rwm <- 0
error_rwm <- 0
final_mix <- 0
final_mix25 <- 0
var_mix <- 0
var_mix25 <- 0
error_mix <- 0
error_mix25 <- 0


e0_true <- 23
emax_true <- 100
e50_true <- 10

o_e0_true <- sqrt(0.5)
o_emax_true <- sqrt(0.3)
o_e50_true <- sqrt(0.3)
a_true = 1
final_mix <- 0
true_param <- data.frame("e0" = e0_true, "emax" = emax_true, "e50" = e50_true, "omega2.e0"=o_e0_true^2 ,"omega2.emax"= o_emax_true^2,"omega2.e50"= o_e50_true^2, "a" =0.1)


seed0 = 39546
replicate = 3


for (j in 1:replicate){

     
# model2 <- inlineModel("
#                       [COVARIATE]
#                       input={p_F}

#                       DEFINITION:
#                       gender = { type        = categorical, 
#                                  categories  = {0,1},
#                                  P(gender=0) = p_F }
                      
#                       [INDIVIDUAL]
#                       input={e0_pop,o_e0,emax_pop,o_emax,gender, beta_F,e50_pop,o_e50}
#                       gender={type=categorical,categories={0,1}}

                      
#                       DEFINITION:
#                       e0  ={distribution=lognormal, 
#                             prediction=e0_pop,
#                             sd=o_e0}
#                       emax   ={distribution=lognormal, prediction=emax_pop,   sd=o_emax}
#                       e50  ={distribution=lognormal, 
#                               reference=e50_pop,  
#                               covariate = gender,
#                               coefficient  = {beta_F,0},
#                               sd=o_e50}

#                        [LONGITUDINAL]
#                       input = {e0, emax, e50,a}

                      
#                       EQUATION:
#                       Cc = e0+emax*t/(e50+t)
                      
#                       DEFINITION:
#                       y1 ={distribution=normal, prediction=Cc, sd=a}

#                       ")


# p <- c(p_F=0.3,  beta_F=0.6,e0_pop=e0_true, o_e0=o_e0_true,
#        emax_pop=emax_true, o_emax=o_emax_true, 
#        e50_pop=e50_true, o_e50=o_e50_true,  
#        a=a_true)
# ind <- list(name=c("gender"))

# adm  <- list(amount=1, time=seq(0,50,by=50))


# y1 <- list(name='y1', time=seq(0,to=50,by=10))


# res2a2 <- simulx(model = model2,
#                  treatment = adm,
#                  parameter = p,
#                  group = list(size=1000, level="covariate"),
#                  output = list(ind, y1))

    # writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_pd.csv")
  # pd.data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_pd.csv", header=T, sep=",")
  # pd.data <- pd.data[pd.data[,4]!=1 ,c(1,2,3,5)]
  # pd.data[,3] <- res2a2$y1[,3]
  # colnames(pd.data) <- c("subject","dose","response","gender")


# saemix.data<-saemixData(name.data=pd.data,header=TRUE,name.group=c("subject"),
# name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
# units=list(x="mg",y="-",covariates="-"))

# saemix.model<-saemixModel(model=modelemax,description="Emax model",type="structural",
# psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,
# dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
# covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
# fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
# byrow=TRUE),error.model="constant")
  
model2 <- inlineModel("
                      
                      [INDIVIDUAL]
                      input={e0_pop,o_e0,emax_pop,o_emax,e50_pop,o_e50}

                      
                      DEFINITION:
                      e0  ={distribution=lognormal, 
                            prediction=e0_pop,
                            sd=o_e0}
                      emax   ={distribution=lognormal, prediction=emax_pop,   sd=o_emax}
                      e50  ={distribution=lognormal, 
                              reference=e50_pop,  
                              sd=o_e50}

                       [LONGITUDINAL]
                      input = {e0, emax, e50,a}

                      
                      EQUATION:
                      Cc = e0+emax*t/(e50+t)
                      
                      DEFINITION:
                      y1 ={distribution=normal, prediction=Cc, sd=a}

                      ")


p <- c(e0_pop=e0_true, o_e0=o_e0_true,
       emax_pop=emax_true, o_emax=o_emax_true, 
       e50_pop=e50_true, o_e50=o_e50_true,  
       a=a_true)

y1 <- list(name='y1', time=seq(0,to=50,by=10))


res2a2 <- simulx(model = model2,
                 parameter = p,
                 group = list(size=1000, level="individual"),
                 output = y1)


  
writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_pd.csv")
  pd.data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_pd.csv", header=T, sep=",")
  pd.data[,3] <- res2a2$y1[,3]
  colnames(pd.data) <- c("subject","time","response")


  saemix.data<-saemixData(name.data=pd.data,header=TRUE,name.group=c("subject"),
name.predictors=c("time"),name.response=c("response"))

saemix.model<-saemixModel(model=modelemax,description="Emax model",type="structural",
psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")
 
  print(j)


  options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100, sampling='')
  theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  ML <- theo_ref[,2:8]
  ML[1:(end+1),]<- theo_ref[end+1,2:8]
  # ML[1:(end+1),1:7]<- true_param
  error_rwm <- error_rwm + (theo_ref[,2:8]-ML)^2
  theo_ref['individual'] <- j
  final_ref <- rbind(final_ref,theo_ref)


  options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50, sampling='seq')
  theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
  theo_mix <- cbind(iterations, theo_mix)

  ML <- theo_mix[,2:8]
  ML[1:(end+1),]<- theo_mix[end+1,2:8]
  # ML[1:(end+1),1:7]<- true_param
  error_mix <- error_mix + (theo_mix[,2:8]-ML)^2
  theo_mix['individual'] <- j
  final_mix <- rbind(final_mix,theo_mix)
  
  options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25, sampling='seq')
  theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
  theo_mix25 <- cbind(iterations, theo_mix25)
  
  ML <- theo_mix25[,2:8]
  ML[1:(end+1),]<- theo_mix25[end+1,2:8]
  # ML[1:(end+1),1:7]<- true_param
  error_mix25 <- error_mix25 + (theo_mix25[,2:8]-ML)^2
  theo_mix25['individual'] <- j
  final_mix25 <- rbind(final_mix25,theo_mix25)
}


graphConvMC_diff(final_ref,final_ref,final_ref)
graphConvMC_diff(final_ref,final_mix,final_mix25)

error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix
error_mix25 <- 1/replicate*error_mix25

err_mix<- theo_ref[-1,]
err_rwm<- theo_ref[-1,]
err_mix25<- theo_ref[-1,]

err_rwm[,2:8] <- error_rwm[-1,1:7]
err_mix[,2:8] <- error_mix[-1,1:7]
err_mix25[,2:8] <- error_mix25[-1,1:7]


err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = seq(1, 4*end, by=4)
err_mix_scaled <- err_mix
err_mix_scaled$iterations = seq(1, 2*end, by=2)
err_mix25$iterations = 1:((K1+K2))


# c <- graphConvMC_se2(err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)])
c <- graphConvMC_sec(err_rwm_scaled[1:end,c(1,3,9)],err_mix_scaled[1:end,c(1,3,9)],err_mix25[1:end,c(1,3,9)])
d <- graphConvMC_sed(err_rwm_scaled[1:end,c(1,6,9)],err_mix_scaled[1:end,c(1,6,9)],err_mix25[1:end,c(1,6,9)])

grid.arrange(c,d, ncol=2)





# c <- graphConvMC_se2(err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)])
e <- graphConvMC_sec(err_rwm_scaled[1:end,c(1,3,9)],err_mix_scaled[1:end,c(1,3,9)],err_mix25[1:end,c(1,3,9)])
f <- graphConvMC_sed(err_rwm_scaled[1:end,c(1,6,9)],err_mix_scaled[1:end,c(1,6,9)],err_mix25[1:end,c(1,6,9)])

grid.arrange(e,f, ncol=2)


