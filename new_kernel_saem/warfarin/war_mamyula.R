# setwd("/Users/karimimohammedbelhal/Desktop/variationalBayes/mcmc_R_isolate/Dir2")
#   source('compute_LL.R') 
#   source('func_aux.R') 
#   source('func_cov.R') 
#   source('func_distcond.R') 
#   source('func_FIM.R') 
#   source('func_ggplot2.R') 
#   source('func_plots.R') 
#   source('func_simulations.R') 
#   source('ggplot2_global.R') 
#   # source('KL.R') 
#   #source('vi.R') 
#   source('global.R')
#   source('main.R')
#   source('mcmc_main.R') 
#   source('main_estep.R')
#   source('main_estep_mcmc.R') 
#   source('main_estep_morekernels.R') 
#   source('main_initialiseMainAlgo.R') 
#   source('main_mstep.R') 
#   source('SaemixData.R')
#   source('plots_ggplot2.R') 
#   source('saemix-package.R') 
#   source('SaemixModel.R') 
#   source('SaemixRes.R') 
#   source('SaemixObject.R') 
#   source('zzz.R') 
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/cmaes/Dir")
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_cov.R') 
  source('func_distcond.R') 
  source('func_FIM.R') 
  source('func_ggplot2.R') 
  source('func_plots.R') 
  source('func_simulations.R') 
  source('ggplot2_global.R') 
  # source('KL.R') 
  #source('vi.R') 
  source('global.R')
  source('main.R')
  source('mcmc_main.R') 
  source('main_estep.R')
  source('main_estep_mcmc.R') 
  source('main_estep_morekernels.R') 
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('plots_ggplot2.R') 
  source('saemix-package.R') 
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem")
source('newkernel_main.R')
source('main_new.R')
source('main_estep_new.R')
source('main_estep_new2.R')
source('main_gd.R')
source('main_estep_gd.R')
source('main_estep_newkernel.R')
source('main_gd_mix.R')
source('main_estep_gd_mix.R')
source('main_estep_mix.R')
source('main_estep_newkernel.R')
source('main_mamyula.R')
source('main_estep_mala.R')
source("mixtureFunctions.R")
library("mlxR")
library(sgd)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")


# Doc
# data(theo.saemix)
# theo.saemix_less <- theo.saemix[1:120,]
# # theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# saemix.data<-saemixData(name.data=theo.saemix_less,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")



setwd("/Users/karimimohammedbelhal/Desktop/CSDA_code/warfarin")
# source('dataproc.R')
# model 
model<-"warfarin_project_model.txt"
# treatment
trt <- read.table("treatment.txt", header = TRUE) 

# parameters 
originalId<- read.table('originalId.txt', header=TRUE) 
populationParameter<- read.vector('populationParameter.txt') 
individualCovariate<- read.table('individualCovariate.txt', header = TRUE) 
list.param <- list(populationParameter,individualCovariate)
# output 
name<-"y1"
time<-read.table("output1.txt",header=TRUE)
out1<-list(name=name,time=time) 
name<-"y2"
time<-read.table("output2.txt",header=TRUE)
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

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem")
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
saemix.model<-saemixModel(model=model1cpt,description="warfarin"
  ,psi0=matrix(c(1,7,1),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE))


saemix.data<-saemixData(name.data=warfarin.saemix_less,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

K1 = 300
K2 = 50
iterations = 1:(K1+K2+1)
gd_step = 0.01
end = K1+K2
seed0 = 395246
#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,0,0,0,0,0), nbiter.saemix = c(K1,K2))
theo_ref<-data.frame(saemix_mamyula(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

theo_ref[end,]

graphConvMC_twokernels(theo_ref,theo_ref, title="new kernel")
#saem with mala
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,2,0,0),nbiter.saemix = c(K1,K2),sigma.val = 0.01,gamma.val=0.01)
theo_mala<-data.frame(saemix_mamyula(saemix.model,saemix.data2,options.mala))
theo_mala <- cbind(iterations, theo_mala)


#saem with mamyula
options.mamyula<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,0,2,0),nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01,lambda.val=0.2)
theo_mamyula<-data.frame(saemix_mamyula(saemix.model,saemix.data,options.mamyula))
theo_mamyula <- cbind(iterations, theo_mamyula)

theo_ref[end,]
theo_mamyula[end,]

graphConvMC_twokernels(theo_ref,theo_mala, title="new kernel")
graphConvMC_threekernels(theo_ref,theo_mala,theo_mamyula, title="new kernel")
graphConvMC_threekernels(theo_ref,theo_mamyula,theo_mamyula, title="new kernel")

#First run on the same dataset

replicate = 3

final_rwm <- 0
final_mix <- 0
for (m in 1:replicate){
  print(m)
  print(m)
  l = list(c(1,5,1,0,0,0),c(0.8,4,0.8,0,0,0),c(1.2,3,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model<-saemixModel(model=model1cpt,description="warfarin"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE))

  options<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0)
  theo_ref<-data.frame(saemix_new_mix(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,theo_ref[-1,])

  options.new<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),map.range=c(1:10),nbiter.burn =0)
  theo_new_ref<-data.frame(saemix_new_mix(saemix.model,saemix.data,options.new))
  theo_mix <- cbind(iterations, theo_new_ref)
  theo_mix['individual'] <- m
  final_mix <- rbind(final_mix,theo_mix[-1,])
}
