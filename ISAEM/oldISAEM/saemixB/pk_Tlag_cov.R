
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
  k<-psi[id,4]
  CL<-k*V
  dt <- pmax(time-Tlag, 0)
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*dt)-exp(-ka*dt))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.2,1,7,1),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","k"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE))



K1 = 200
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50
#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
warfa<-cbind(iterations,warfa)

options_warfaincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='seq')
warfaincr25<-data.frame(saemix_incremental(saemix.model_warfa,saemix.data_warfa,options_warfaincr25))
warfaincr25<-cbind(iterations,warfaincr25)

options_warfaincr50<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
warfaincr50<-data.frame(saemix_incremental(saemix.model_warfa,saemix.data_warfa,options_warfaincr50))
warfaincr50<-cbind(iterations,warfaincr50)

graphConvMC2_saem(warfa,warfa, title="new kernel") 

warfa$algo <- 'rwm'
warfaincr25$algo <- 'ISAEM25'
warfaincr50$algo <- 'ISAEM50'

warfa_scaled <- warfa[rep(seq_len(nrow(warfa)), each=100/batchsize25),]
warfa_scaled$iterations = 1:(4*(K1+K2+1))

warfa_scaled50 <- warfaincr50[rep(seq_len(nrow(warfaincr50)), each=100/batchsize50),]
warfa_scaled50$iterations = 1:(2*(K1+K2+1))

comparison <- 0
# comparison <- rbind(theo_ref,theo_incremental)
# comparison <- rbind(warfa_scaled[iterations,],warfaincr)
comparison <- rbind(warfa_scaled[iterations,],warfa_scaled50[iterations,],warfaincr25)

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(1:3))
warfanew<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew))


#Diff initial values

replicate = 3
seed0 = 395246

#RWM
final_rwm <- 0
final_50 <- 0
final_25 <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(1,7,1,0,0,0),c(0.8,7.2,0.8,0,0,0),c(1.2,6.8,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

  options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
  warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
  warfa<-cbind(iterations,warfa)
  warfa['individual'] <- m
  warfa$algo <- 'rwm'
  warfa_scaled <- warfa[-1,]
  warfa_scaled$iterations = seq(1, 4*end, by=4)
  warfa_scaled <- warfa_scaled[rep(seq_len(nrow(warfa_scaled)), each=100/batchsize25),]
  final_rwm <- rbind(final_rwm,warfa_scaled[0:end,])

  options_warfaincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
  warfaincr25<-data.frame(saemix_incremental(saemix.model_warfa,saemix.data_warfa,options_warfaincr25))
  warfaincr25<-cbind(iterations,warfaincr25)
  warfaincr25['individual'] <- m
  warfaincr25$algo <- 'ISAEM25'
  warfaincr25 <- warfaincr25[-1,]
  warfaincr25$iterations = seq(1, end, by=1)
  final_25 <- rbind(final_25,warfaincr25[0:end,])

  options_warfaincr50<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
  warfaincr50<-data.frame(saemix_incremental(saemix.model_warfa,saemix.data_warfa,options_warfaincr50))
  warfaincr50<-cbind(iterations,warfaincr50)
  warfaincr50['individual'] <- m
  warfaincr50$algo <- 'ISAEM50'
  warfa_scaled50 <- warfaincr50[-1,]
  warfa_scaled50$iterations = seq(1, 2*end, by=2)
  warfa_scaled50 <- warfa_scaled50[rep(seq_len(nrow(warfa_scaled50)), each=100/batchsize50),]
  final_50 <- rbind(final_50,warfa_scaled50[0:end,])
}

a <- graphConvMC_diffz(final_rwm[,c(1,3,9)],final_50[,c(1,3,9)],final_25[,c(1,3,9)])
b <- graphConvMC_diffw(final_rwm[,c(1,6,9)],final_50[,c(1,6,9)],final_25[,c(1,6,9)])

grid.arrange(a,b, ncol=2)




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
k_true <- 0.1/8
o_Tlag <- 0.57
o_ka <- 0.5
o_V <- 0.2
o_k <- 0.3

true_param <- data.frame("Tlag" = Tlag_true,"ka" = ka_true, "V" = V_true, "k" = k_true, "omega2.Tlag"=o_Tlag^2, "omega2.ka"=o_ka^2 ,"omega2.V"= o_V^2,"omega2.k"= o_k^2, "a" = 1)
seed0 = 39546
replicate = 20
for (m in 1:replicate){
  
    model<-"/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/warfarin_project_model.txt"
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
name<-"y2"
time<-read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/output2.txt",header=TRUE)
out2<-list(name=name,time=time) 
out<-list(out1,out2)

# call the simulator 
res <- simulx(model=model,treatment=trt,parameter=list.param,output=out)
# warfarin.saemix <- data(res)

# writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/warfarin/war_synth.csv")
# table <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/warfarin/war_synth.csv", header=T, sep=",")
# table <- table[table$ytype==1,]

# table[,5] <- 0

warfarin.saemix <- res$y1
warfarin.saemix["amount"] <- 0
treat <- res$treatment
treat["y1"] <- 0
treat <- treat[c(1,2,4,3)]

j <- 1
l<-c()


for (i in 1:241) {
    
    if(t(warfarin.saemix["id"])[i]==t(treat["id"])[j]){
        print(rownames(warfarin.saemix[i,]))
        l <- rbind(l,rownames(warfarin.saemix[i,]))
        j<-j+1
      }
}

warfarin.saemix <- rbind(treat[1,], warfarin.saemix)
warfarin.saemix[1:7,4] <- treat[1,4]
j <- 2
for (i in l[-1]){
  warfarin.saemix[(as.numeric(i)+1):(as.numeric(i)+length(which(t(warfarin.saemix["id"]==j)))),4] <- treat[j,4]
  warfarin.saemix <- rbind(warfarin.saemix[1:(as.numeric(i)-1),], treat[j,], warfarin.saemix[(as.numeric(i)+1):nrow(warfarin.saemix),])
  j <- j +1
}

rownames(warfarin.saemix) <- 1:nrow(warfarin.saemix)
# warfarin.saemix <- table[c(1,2,3,5)]

saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.2,3,13,2),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","k"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE))

saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")



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

#RANDOMPASS
options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randompass')
theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
theo_mix <- cbind(iterations, theo_mix)
ML <- theo_mix[,2:10]
# ML[1:(end+1),]<- theo_mix[end+1,2:10]
ML[1:(end+1),1:9]<- true_param
error_mixpass <- error_mixpass + (theo_mix[,2:10]-ML)^2
theo_mix['individual'] <- m
final_mix <- rbind(final_mix,theo_mix)

options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randompass')
theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))  
theo_mix25 <- cbind(iterations, theo_mix25)
ML <- theo_mix25[,2:10]
# ML[1:(end+1),]<- theo_mix25[end+1,2:10]
ML[1:(end+1),1:9]<- true_param
error_mix25pass <- error_mix25pass + (theo_mix25[,2:10]-ML)^2
theo_mix25['individual'] <- m
final_mix25 <- rbind(final_mix25,theo_mix25)

#RANDOMITER
options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50,sampling='randomiter')
theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
theo_mix <- cbind(iterations, theo_mix)
ML <- theo_mix[,2:10]
# ML[1:(end+1),]<- theo_mix[end+1,2:10]
ML[1:(end+1),1:9]<- true_param
error_mixiter <- error_mixiter + (theo_mix[,2:10]-ML)^2
theo_mix['individual'] <- m
final_mix <- rbind(final_mix,theo_mix)

options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25,sampling='randomiter')
theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_mix25)
ML <- theo_mix25[,2:10]
# ML[1:(end+1),]<- theo_mix25[end+1,2:10]
ML[1:(end+1),1:9]<- true_param
error_mix25iter <- error_mix25iter + (theo_mix25[,2:10]-ML)^2
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

err_rwm[,2:10] <- error_rwm[-1,]

err_mixseq[,2:10] <- error_mixseq[-1,]
err_mix25seq[,2:10] <- error_mix25seq[-1,]

err_mixpass[,2:10] <- error_mixpass[-1,]
err_mix25pass[,2:10] <- error_mix25pass[-1,]

err_mixiter[,2:10] <- error_mixiter[-1,]
err_mix25iter[,2:10] <- error_mix25iter[-1,]


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
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,12,13)],err_mixseq_scaled [0:end,c(1,i,12,13)],err_mix25seq[0:end,c(1,i,12,13)],
                                              err_mixpass_scaled [0:end,c(1,i,12,13)],err_mix25pass[0:end,c(1,i,12,13)],
                                              err_mixiter_scaled [0:end,c(1,i,12,13)],err_mix25iter[0:end,c(1,i,12,13)])

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


