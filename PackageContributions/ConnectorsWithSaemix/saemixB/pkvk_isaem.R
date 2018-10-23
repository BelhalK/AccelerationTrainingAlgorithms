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
  tauD<-psi[id,1]
  tauS<-psi[id,2]
  ypred<-1
  return(ypred)
}

pkvk_data <- readDatamlx(project = project.file)
treat <- pkvk_data$treatment[,c(3)]
# pkvk.saemix <- cbind(pkvk_data$Y,treat)
pkvk.saemix <- pkvk_data$y

saemix.data<-saemixData(name.data=pkvk.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("time"),name.response=c("y"), name.X="time")

saemix.model<-saemixModel(model=model1cpt,description="pkvk",type="structural"
  ,psi0=matrix(c(10,10),ncol=2,byrow=TRUE, dimnames=list(NULL, c("tauD","tauS"))),
  fixed.estim=c(1,1),transform.par=c(1,1),omega.init=matrix(diag(2),ncol=2,byrow=TRUE),
  covariance.model=matrix(diag(2),ncol=2,byrow=TRUE), error.model = "proportional")

K1 = 1000
K2 = 500
iterations = 1:(K1+K2+1)
end = K1+K2

runtime = 20

options<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=1,monolix=TRUE,
 nb.replacement=100,sampling='randompass', duration = runtime)
pkvk<-data.frame(saemix(saemix.model,saemix.data,options))
pkvk <- cbind(iterations, pkvk)
row_sub_ref  = apply(pkvk, 1, function(row) all(row !=0 ))
pkvk <- pkvk[row_sub_ref,]
pkvk$algo <- 'full'
pkvk$iterations <- seq(0,runtime, length.out=length(pkvk$iterations))



options50<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=1,monolix=TRUE,
 nb.replacement=50,sampling='randompass', duration = runtime)
pkvk50<-data.frame(saemix(saemix.model,saemix.data,options50))
pkvk50 <- cbind(iterations, pkvk50)
row_sub_50  = apply(pkvk50, 1, function(row) all(row !=0 ))
pkvk50 <- pkvk50[row_sub_50,]
pkvk50$algo <- 'half'
pkvk50$iterations <- seq(0,runtime, length.out=length(pkvk50$iterations))


options25<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=1,monolix=TRUE,
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
# comparison <- rbind(pkvk,pkvk50)
comparison <- rbind(pkvk,pkvk25,pkvk50)
var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
prec <- seplot(var, title="comparison",legend=TRUE)


