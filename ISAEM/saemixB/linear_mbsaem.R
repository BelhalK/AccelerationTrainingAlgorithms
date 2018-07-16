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
# load("isaem_linear.RData")
# save.image("isaem_linear.RData")

library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

K1 = 400
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

# Doc
oxboys.saemix<-read.table( "data/ox_synth.csv",header=T,na=".",sep=",")
oxboys.saemix_less <- oxboys.saemix[,]
n <- length(unique(oxboys.saemix_less$id))

saemix.data<-saemixData(name.data=oxboys.saemix_less,header=TRUE,
  name.group=c("id"),name.predictors=c("time"),name.response=c("y"),
  units=list(x="yr",y="cm"))

# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

growth.linear<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (2 columns, base and slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  base<-psi[id,1]
  slope<-psi[id,2]
  f<-base+slope*x
  return(f)
}

saemix.model<-saemixModel(model=growth.linear,description="Linear model",type="structural",
  psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(0,0),covariance.model=matrix(c(1,0,0,1),ncol=2,byrow=TRUE),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE), 
  error.model="constant")



K1 = 300
K2 = 100
iterations = 0:(K1+K2-1)
end = K1+K2
batchsize25 = 25
batchsize50 = 50

seed0=3456

options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100,sampling='seq')
theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options)$param)
theo_ref <- cbind(iterations, theo_ref[-1,])


options.incremental50<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(0,0,0,6), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,
                          nbiter.burn =0, nb.replacement=50,sampling='randompass')
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
