setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incrementalwithMlxConnectors/R")
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


  source('/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incrementalwithMlxConnectors/R/mixtureFunctions.R')
  source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incrementalwithMlxConnectors/plots.R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/incrementalwithMlxConnectors")

library("rlist")
# library("mlxR")
library(MlxConnectors)
initializeMlxConnectors(software = "monolix")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
require(madness)


project.file <- "mlxproject_ode/test_ode.mlxtran"
loadProject(project.file)
getEstimatedPopulationParameters()


computePredictions(getEstimatedIndividualParameters()$saem)
computePredictions(getEstimatedIndividualParameters()$saem, individualIds = c(10,20))

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  time<-xidep[,2]  
  r<-psi[id,1]
  alpha<-psi[id,2]
  V<-psi[id,3]
  Cl<-psi[id,4]
  k <- Cl/V
  ypred<-r + alpha + V + Cl
  return(ypred)
}



library("mlxR")



warfa_data <- readDatamlx(project = project.file)
treat <- warfa_data$treatment[,c(3)]
warfa.saemix <- cbind(warfa_data$y,treat)
saemix.data<-saemixData(name.data=warfa.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("treat","time"),name.response=c("y"), name.X="time")


# warfa_data <- read.table("mlxproject_ode/oral_data_saemix.txt", header=T)
# saemix.data<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
#   name.predictors=c("treat","time"),name.response=c("y"), name.X="time")

# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,1,1,1),ncol=4,byrow=TRUE, dimnames=list(NULL, c("r","alpha","V","Cl"))),fixed.estim=c(1,1,1,1),
  transform.par=c(1,1,1,1),omega.init=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),ncol=4,byrow=TRUE))


K1 = 100
K2 = 50
iterations = 1:(K1+K2)
end = K1+K2
batchsize25 = 25
batchsize50 = 50

seed0=3456


options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1,nbiter.mcmc = c(2,2,2,0),
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, 
 map.range=c(0), nb.replacement=100,sampling='seq')
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])
row_sub_ref  = apply(theo_ref, 1, function(row) all(row !=0 ))
theo_ref <- theo_ref[row_sub_ref,]
theo_ref$algo <- 'full'
theo_ref$iterations <- seq(0,10, length.out=length(theo_ref$iterations))



options<-list(seed=39546,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1,nbiter.mcmc = c(2,2,2,2),
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, 
 map.range=c(1:4), nb.replacement=100,sampling='seq')
theo_ref<-saemix_incremental(saemix.model,saemix.data,options)
theo_ref <- data.frame(theo_ref$param)
theo_ref <- cbind(iterations, theo_ref[-1,])
row_sub_ref  = apply(theo_ref, 1, function(row) all(row !=0 ))
theo_ref <- theo_ref[row_sub_ref,]
theo_ref$algo <- 'full'
theo_ref$iterations <- seq(0,10, length.out=length(theo_ref$iterations))


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(2,2,2,2), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(1:4),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo_mix25<-saemix_incremental(saemix.model,saemix.data,options.incremental25)
theo_mix25 <- data.frame(theo_mix25$param)
theo_mix25 <- cbind(iterations, theo_mix25[-1,])
row_sub  = apply(theo_mix25, 1, function(row) all(row !=0 ))
theo_mix25 <- theo_mix25[row_sub,]
theo_mix25$algo <- 'quarter'
theo_mix25$iterations <- seq(0,10, length.out=length(theo_mix25$iterations))


comparison <- 0
comparison <- rbind(theo_ref[,],theo_mix25[,])
var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)


seplot <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  graf <- ggplot(df)+geom_line(aes(iterations,value,by=value,colour = df$algo),show.legend = legend) +
  xlab("iterations") + ylab('value') + facet_wrap(~variable,scales = "free_y") + theme_bw() 
  grid.arrange(graf)
  # do.call("grid.arrange", c(graf, ncol=1, top=title))
}


prec <- seplot(var, title="comparison",legend=TRUE)
assign(paste("prec", i, sep = ""), prec) 

