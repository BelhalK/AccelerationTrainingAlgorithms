# savedwarfarin.saemix <- warfarin.saemix
# load("pkcov_samplingstrat.RData")
# save.image("pkcov_samplingstrat.RData")
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
# load("isaem_design2.RData")
# save.image("isaem_design2.RData")

library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)



# warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/warfarin_data.txt", header=T)
# saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")



# saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("amount","time"),name.response=c("y1"), name.X="time", name.covariates=c("wt"),units=list(x="kg",
#   covariates=c("kg/ha")))

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



saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time", name.covariates=c("wt"),units=list(x="kg",
  covariates=c("kg/ha")))



saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.2,3,10,2),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","Cl"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE),fixed.estim=c(1,0,0,0),covariate.model=t(c(0,0,1,1)),error.model="constant")


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


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass', gamma=1)
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])

summary <- theo50$summary
chosen <- data.frame(theo50$chosen)
kiter <- 30
test <- t(summary)
test <- data.frame(test)
test$iterations <- 1:kiter
df <- melt(test ,  id.vars = 'iterations')
current <- theo_mix50[1:kiter,2]
chosen <- t(chosen)
chosen <- data.frame(chosen)
chosen$iterations <- 1:kiter
df.chosen <- melt(chosen ,  id.vars = 'iterations')
df$chosen <- df.chosen$value

ggplot(df, aes(iterations,value)) + geom_point(aes(colour = chosen))+ 
  geom_point(data = theo_mix50[1:kiter,], aes(x = iterations, y = theo_mix50[end,2]), color = "red")+ 
  geom_point(data = theo_mix50[1:kiter,], aes(x = iterations, y = theo_mix50[1:kiter,2]), color = "yellow")+ theme_bw()



options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass', gamma=0)
theo50_random<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50_random <- data.frame(theo50_random$param)
theo_mix50_random <- cbind(iterations, theo_mix50_random[-1,])

summary_random <- theo50_random$summary
chosen_random <- data.frame(theo50_random$chosen)
kiter <- 30
test <- t(summary_random)
test <- data.frame(test)
test$iterations <- 1:kiter
df <- melt(test ,  id.vars = 'iterations')
current <- theo_mix50_random[1:kiter,2]
chosen_random <- t(chosen_random)
chosen_random <- data.frame(chosen_random)
chosen_random$iterations <- 1:kiter
df.chosen_random <- melt(chosen_random ,  id.vars = 'iterations')
df$chosen_random <- df.chosen_random$value


ggplot(df, aes(iterations,value)) + geom_point(aes(colour = chosen_random))+ 
  geom_point(data = theo_mix50_random[1:kiter,], aes(x = iterations, y = theo_mix50_random[end,2]), color = "red")+ 
  geom_point(data = theo_mix50_random[1:kiter,], aes(x = iterations, y = theo_mix50_random[1:kiter,2]), color = "yellow")+ theme_bw()


theo_ref_scaled <- theo_ref
theo_mix50_scaled <- theo_mix50
theo_mix50_random_scaled <- theo_mix50_random


theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
theo_mix50_random_scaled$iterations = theo_mix50_random_scaled$iterations*0.5

graphConvMC_threekernels(theo_ref_scaled,theo_mix50_scaled,theo_mix50_random_scaled)





options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
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

graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix50_scaled,theo_mix50_scaled)



###NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL####
options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(1:5), nb.replacement=100,sampling='seq')
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo50<-saemix_incremental(saemix.model,saemix.data,options.incremental50)
theo_mix50 <- data.frame(theo50$param)
theo_mix50 <- cbind(iterations, theo_mix50[-1,])

options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
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

graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix50_scaled,theo_mix50_scaled)

graphConvMC_threekernels(theo_mix25_scaled,theo_mix50_scaled,theo_mix25_scaled)
###NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL#######NEW KERNEL####


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

graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix50_scaled,theo_mix75_scaled,theo_mix85_scaled)

K1 = 600
K2 = 100
iterations = 0:(K1+K2-1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50


final_ref <- 0
final_mix50 <- 0
final_mix25 <- 0
final_mix75 <- 0
final_mix85 <- 0

error_rwm <- 0
error_mix50 <- 0
error_mix25 <- 0
error_mix75 <- 0
error_mix85 <- 0


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
replicate = 40
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

trt <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/design2/treatment.txt", header = TRUE) 
originalId<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/design2/originalId.txt', header=TRUE) 
individualCovariate<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/design2/individualCovariate.txt', header = TRUE) 
time<-read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/design2/output1.txt",header=TRUE)

# trt <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/treatment.txt", header = TRUE) 
# originalId<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/originalId.txt', header=TRUE) 
# individualCovariate<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/individualCovariate.txt', header = TRUE) 
# time<-read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/output1.txt",header=TRUE)


list.param <- list(populationParameter,individualCovariate)
name<-"y1"
out1<-list(name=name,time=time) 

# call the simulator 
res <- simulx(model=myModel,treatment=trt,parameter=list.param,output=out1)
# writeDatamlx(res, result.file = paste("/Users/karimimohammedbelhal/Desktop/data_pk/pk_mcstudy_", m, ".csv", sep=""))
# writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Desktop/data_pk/isaem_design2.csv")
# warfarin.saemix<-read.table(paste("/Users/karimimohammedbelhal/Desktop/data_pk/pk_mcstudy_", m, ".csv", sep=""), header=T, sep=",")
# typeof(warfarin.saemix[2,3])
# typeof(res$y1[2,3])
individualCovariate$wt <- log(individualCovariate$wt/70)
warfarin.saemix <- res$y1
treat <- res$treatment[,c(1,3)]
covandtreat <- merge(individualCovariate ,treat,by="id")
warfarin.saemix <- merge(covandtreat ,warfarin.saemix,by="id")

 
# warfarin.saemix <- table[c(1,2,3,5)]

# saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(0.2,3,10,2),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","Cl"))),
#   transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
#   byrow=TRUE),error.model="constant")

# saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.2,3,10,2),ncol=4,byrow=TRUE, dimnames=list(NULL, c("Tlag","ka","V","Cl"))),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4, 
  byrow=TRUE),covariate.model=matrix(c(0,0,1,1),ncol=4,byrow=TRUE),error.model="constant")

saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time", name.covariates=c("wt"),units=list(x="kg",
  covariates=c("kg/ha")))



options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='')
theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref[-1,])

ML <- theo_ref[,2:10]
# ML[1:(end+1),]<- theo_ref[end+1,2:10]
ML[1:nrow(ML),]<-true_param

error_rwm <- error_rwm + (theo_ref[,2:10]-ML)^2
theo_ref['individual'] <- m
final_ref <- rbind(final_ref,theo_ref)



options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo_mix50<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental50))
theo_mix50 <- cbind(iterations, theo_mix50[-1,])

# ML[1:(end+1),]<- theo_mix50[end+1,2:10]

error_mix50 <- error_mix50 + (theo_mix50[,2:10]-ML)^2
theo_mix50['individual'] <- m
final_mix50 <- rbind(final_mix50,theo_mix50)


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_mix25[-1,])

# ML[1:(end+1),]<- theo_mix25[end+1,2:10]

error_mix25 <- error_mix25 + (theo_mix25[,2:10]-ML)^2
theo_mix25['individual'] <- m
final_mix25 <- rbind(final_mix25,theo_mix25)

options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
theo_mix75<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental75))
theo_mix75 <- cbind(iterations, theo_mix75[-1,])

# ML[1:(end+1),]<- theo_mix75[end+1,2:10]

error_mix75 <- error_mix75 + (theo_mix75[,2:10]-ML)^2
theo_mix75['individual'] <- m
final_mix75 <- rbind(final_mix75,theo_mix75)

options.incremental85<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=85,sampling='randompass')
theo_mix85<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental85))
theo_mix85 <- cbind(iterations, theo_mix85[-1,])

# ML[1:(end+1),]<- theo_mix85[end+1,2:10]

error_mix85 <- error_mix85 + (theo_mix85[,2:10]-ML)^2
theo_mix85['individual'] <- m
final_mix85 <- rbind(final_mix85,theo_mix85)
 
}


graphConvMC_diff(final_ref,final_mix75,final_mix85)



final_ref$algo <- '100'
final_mix50$algo <- '50'
final_mix25$algo <- '25'
final_mix75$algo <- '75'
final_mix85$algo <- '85'


final_ref_scaled <- final_ref
final_mix50_scaled <- final_mix50
final_mix25_scaled <- final_mix25
final_mix75_scaled <- final_mix75
final_mix85_scaled <- final_mix85


final_ref_scaled$iterations = final_ref_scaled$iterations*1
final_mix50_scaled$iterations = final_mix50_scaled$iterations*0.5
final_mix25_scaled$iterations = final_mix25_scaled$iterations*0.25
final_mix75_scaled$iterations = final_mix75_scaled$iterations*0.75
final_mix85_scaled$iterations = final_mix85_scaled$iterations*0.85



for (i in c(3,7)){
comparison <- 0

comparison <- rbind(final_ref_scaled[,c(1,i,11,12)],final_mix25_scaled[,c(1,i,11,12)],
                final_mix50_scaled[,c(1,i,11,12)],final_mix75_scaled[,c(1,i,11,12)],
                final_mix85_scaled[,c(1,i,11,12)])

var <- melt(comparison, id.var = c('iterations','algo','individual'), na.rm = TRUE)

beta0 <- seplot(var,colnames(final_ref_scaled)[i], title="comparison",legend=TRUE)

}


for (i in c(3,7)){
comparison <- 0

comparison <- rbind(final_ref_scaled[,c(1,i,11,12)],
                final_mix85_scaled[,c(1,i,11,12)])

var <- melt(comparison, id.var = c('iterations','algo','individual'), na.rm = TRUE)

beta0 <- seplot(var,colnames(final_ref_scaled)[i], title="comparison",legend=TRUE)

}



error_rwm <- 1/replicate*error_rwm
err_rwm<- theo_ref[,]
err_rwm[,2:10] <- error_rwm[,]
err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = err_rwm_scaled$iterations*1
err_rwm_scaled$algo <- 'SAEM'
err_rwm_scaled$method <- 'seq'


error_mix50 <- 1/replicate*error_mix50
err_mix50<- theo_ref[,]
err_mix50[,2:10] <- error_mix50[,]
err_mix50_scaled <- err_mix50
err_mix50_scaled$iterations = err_mix50_scaled$iterations*0.5
err_mix50_scaled$algo <- 'ISAEM 50'
err_mix50_scaled$method <- 'randompass'

error_mix25 <- 1/replicate*error_mix25
err_mix25<- theo_ref[,]
err_mix25[,2:10] <- error_mix25[,]
err_mix25_scaled <- err_mix25
err_mix25_scaled$iterations = err_mix25_scaled$iterations*0.25
err_mix25_scaled$algo <- 'ISAEM 25'
err_mix25_scaled$method <- 'randompass'

error_mix75 <- 1/replicate*error_mix75
err_mix75<- theo_ref[,]
err_mix75[,2:10] <- error_mix75[,]
err_mix75_scaled <- err_mix75
err_mix75_scaled$iterations = err_mix75_scaled$iterations*0.75
err_mix75_scaled$algo <- 'ISAEM 75'
err_mix75_scaled$method <- 'randompass'

error_mix85 <- 1/replicate*error_mix85
err_mix85<- theo_ref[,]
err_mix85[,2:10] <- error_mix85[,]
err_mix85_scaled <- err_mix85
err_mix85_scaled$iterations = err_mix85_scaled$iterations*0.85
err_mix85_scaled$algo <- 'ISAEM 85'
err_mix85_scaled$method <- 'randompass'





for (i in 2:10){
# i = 6
comparison <- 0
comparison <- rbind(err_rwm_scaled[0:end,c(1,i,12,13)],
                    err_mix50_scaled[0:end,c(1,i,12,13)],
                    err_mix25_scaled[0:end,c(1,i,12,13)],
                    err_mix75_scaled[0:end,c(1,i,12,13)],
                    err_mix85_scaled[0:end,c(1,i,12,13)])

var <- melt(comparison, id.var = c('iterations','algo','method'), na.rm = TRUE)


prec <- seplot(var, title="ALGO - EM (same complexity)",legend=TRUE)
# setwd("/Users/karimimohammedbelhal/Desktop/")
# ggsave(paste("precwarfa_", i, ".png", sep=""),prec)
}

