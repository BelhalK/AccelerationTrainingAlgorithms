
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
  source("/Users/karimimohammedbelhal/Desktop/papers/iem_code/imcem_saemix/plots_se.R")
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
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]
  CL<-k*V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))



K1 = 300
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50
#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
warfa<-cbind(iterations,warfa)

options_warfaincr25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
warfaincr25<-data.frame(saemix_incremental(saemix.model_warfa,saemix.data_warfa,options_warfaincr25))
warfaincr25<-cbind(iterations,warfaincr25)

options_warfaincr50<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
warfaincr50<-data.frame(saemix_incremental(saemix.model_warfa,saemix.data_warfa,options_warfaincr50))
warfaincr50<-cbind(iterations,warfaincr50)

graphConvMC2_saem(warfa,warfaincr, title="new kernel")

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


ka_true <- 1
V_true <- 8
k_true <- 0.1/8
o_ka <- 0.5
o_V <- 0.2
o_k <- 0.3

true_param <- c(ka_true,V_true,k_true,o_ka,o_V,o_k)

seed0 = 39546
replicate = 20
replicate <- 20
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
warfarin.saemix_less <- warfarin.saemix[,]
  saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


saemix.data<-saemixData(name.data=warfarin.saemix_less,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")



options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100)
theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)
var_rwm <- var_rwm + (theo_ref[,2:8]-true_param)^2
ML <- theo_ref[,2:8]
ML[1:(end+1),]<- theo_ref[end+1,2:8]
error_rwm <- error_rwm + (theo_ref[,2:8]-ML)^2
theo_ref['individual'] <- j
final_ref <- rbind(final_ref,theo_ref)


options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
theo_mix <- cbind(iterations, theo_mix)

var_mix <- var_mix + (theo_mix[,2:8]-true_param)^2
ML <- theo_mix[,2:8]
ML[1:(end+1),]<- theo_mix[end+1,2:8]
error_mix <- error_mix + (theo_mix[,2:8]-ML)^2
theo_mix['individual'] <- j
final_mix <- rbind(final_mix,theo_mix)

options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_mix25)

var_mix25 <- var_mix25 + (theo_mix25[,2:8]-true_param)^2
ML <- theo_mix25[,2:8]
ML[1:(end+1),]<- theo_mix25[end+1,2:8]
error_mix25 <- error_mix25 + (theo_mix25[,2:8]-ML)^2
theo_mix25['individual'] <- j
final_mix25 <- rbind(final_mix25,theo_mix25)
 
}

graphConvMC_diff(final_ref,final_ref,final_ref)
graphConvMC_diff(final_ref,final_mix,final_mix25)

error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix
error_mix25 <- 1/replicate*error_mix25

error_rwm <- cbind(iterations, error_rwm)
error_mix <- cbind(iterations, error_mix)
error_mix25 <- cbind(iterations, error_mix25)

err_mix<- theo_ref
err_rwm<- theo_ref
err_mix25<- theo_ref


err_rwm[,2:8] <- error_rwm[,2:8]
err_mix[,2:8] <- error_mix[,2:8]
err_mix25[,2:8] <- error_mix25[,2:8]

err_mix[2,] = err_rwm[2,]=err_mix25[2,]

err_rwm_scaled <- err_rwm[rep(seq_len(nrow(err_rwm)), each=4),]
err_mix_scaled <- err_mix[rep(seq_len(nrow(err_mix)), each=2),]
err_rwm_scaled$iterations = 1:(4*(K1+K2+1))
err_mix_scaled$iterations = 1:(2*(K1+K2+1))


# c <- graphConvMC_se2(err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)])
c <- graphConvMC_sec(err_rwm_scaled[2:end,c(1,2,9)],err_mix_scaled[2:end,c(1,2,9)],err_mix25[2:end,c(1,2,9)])
d <- graphConvMC_sed(err_rwm_scaled[2:end,c(1,5,9)],err_mix_scaled[2:end,c(1,5,9)],err_mix25[2:end,c(1,5,9)])

grid.arrange(c,d, ncol=2)
