# setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/PackageContributions/ConnectorsWithSaemix/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/PackageContributions/ConnectorsWithSaemix/saemixB/incrementalR")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('main.R')
  source('main_estep.R')
  source('main_estep_incremental.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/PackageContributions/ConnectorsWithSaemix/saemixB")
library("rlist")
library("mlxR")
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
################################################################ MonolixProject ####################################################################################################################################

project.file <- "mlxProjects/pkvk/pkvk_project.mlxtran"
loadProject(project.file)

# getEstimatedPopulationParameters()
# getEstimatedLogLikelihood()
# runLogLikelihoodEstimation(linearization = FALSE, wait = TRUE)
# computePredictions(getEstimatedIndividualParameters()$saem)
# computePredictions(getEstimatedIndividualParameters()$saem, individualIds = c(10,20))

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  Tlag<-psi[id,1]
  Tk0<-psi[id,2]
  Vol<-psi[id,3]
  Cl<-psi[id,4]
  ke0<-psi[id,5]
  IC50<-psi[id,6]
  gamma<-psi[id,7]
  s<-psi[id,13]
  d<-psi[id,9]
  beta<-psi[id,10]
  delta<-psi[id,11]
  p<-psi[id,12]
  c<-psi[id,13]
  ypred<-1
  return(ypred)
}

pkvk_data <- readDatamlx(project = project.file)
treat <- pkvk_data$treatment[,c(3)]
# pkvk.saemix <- cbind(pkvk_data$Y,treat)
pkvk.saemix <- pkvk_data$y_2

saemix.data<-saemixData(name.data=pkvk.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("time"),name.response=c("y_2"), name.X="time")

# saemix.model<-saemixModel(model=model1cpt,description="pkvk",type="structural"
#   ,psi0=matrix(c(1000,1,0.00005,0.05,20,5,0.9,0.7),ncol=13,byrow=TRUE, dimnames=list(NULL, c("s","d","beta","delta","p","c","eta","epsilon"))),fixed.estim=c(1,1,1,1,1,1,1,1),
#   transform.par=c(1,1,1,1,1,1,3,3),omega.init=matrix(diag(13),ncol=13,byrow=TRUE),covariance.model=matrix(diag(13),ncol=13,byrow=TRUE))

cov.model = matrix(0,nrow=13,ncol=13,byrow=TRUE)
cov.model[1,1] <- 1
cov.model[2,2] <- 1

saemix.model<-saemixModel(model=model1cpt,description="pkvk",type="structural"
  ,psi0=matrix(c(0.01,0.2,1,15,1,1.5,2,10000,0.5,0.00005,1,20,2),ncol=13,byrow=TRUE,
   dimnames=list(NULL, c("Tlag","Tk0","Vol","Cl","ke0","IC50","gamma","s","d","beta","delta","p","c"))),
  fixed.estim=c(1,1,0,0,0,0,0,0,0,0,0,0,0),
  transform.par=c(1,1,1,1,1,1,1,1,1,1,1,1,1),omega.init=matrix(diag(13),ncol=13,byrow=TRUE),
  covariance.model=cov.model)


K1 = 2000
K2 = 1000
iterations = 1:(K1+K2+1)
end = K1+K2

runtime = 180
nchains=3

options<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=100,sampling='randompass', duration = runtime)
pkvk<-data.frame(saemix(saemix.model,saemix.data,options))
pkvk <- cbind(iterations, pkvk)
row_sub_ref  = apply(pkvk, 1, function(row) all(row !=0 ))
pkvk <- pkvk[row_sub_ref,]
pkvk$algo <- 'full'
pkvk$iterations <- seq(0,runtime, length.out=length(pkvk$iterations))



# options75<-list(seed=39546,map=F,fim=F,ll.is=F,
#   nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
#   displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
#  nb.replacement=75,sampling='randompass', duration = runtime)
# pkvk75<-data.frame(saemix(saemix.model,saemix.data,options75))
# pkvk75 <- cbind(iterations, pkvk75)
# row_sub_ref  = apply(pkvk75, 1, function(row) all(row !=0 ))
# pkvk75 <- pkvk75[row_sub_ref,]
# pkvk75$algo <- '75'
# pkvk75$iterations <- seq(0,runtime, length.out=length(pkvk75$iterations))


options50<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=50,sampling='randompass', duration = runtime)
pkvk50<-data.frame(saemix(saemix.model,saemix.data,options50))
pkvk50 <- cbind(iterations, pkvk50)
row_sub_50  = apply(pkvk50, 1, function(row) all(row !=0 ))
pkvk50 <- pkvk50[row_sub_50,]
pkvk50$algo <- 'half'
pkvk50$iterations <- seq(0,runtime, length.out=length(pkvk50$iterations))


options25<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=25,sampling='randompass', duration = runtime)
pkvk25<-data.frame(saemix(saemix.model,saemix.data,options25))
pkvk25 <- cbind(iterations, pkvk25)
row_sub_25  = apply(pkvk25, 1, function(row) all(row !=0 ))
pkvk25 <- pkvk25[row_sub_25,]
pkvk25$algo <- 'quarter'
pkvk25$iterations <- seq(0,runtime, length.out=length(pkvk25$iterations))


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


comparison <- 0
comparison <- rbind(pkvk,pkvk)
# comparison <- rbind(pkvk,pkvk25,pkvk50)
# comparison <- rbind(pkvk,pkvk25,pkvk50,pkvk75)
var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
prec <- seplot(var, title="comparison",legend=TRUE)


