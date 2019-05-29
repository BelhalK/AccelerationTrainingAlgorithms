# save.image("arthritis_10chains.RData")
# save.image("arthritis_40chains.RData")
# load("arthritis_40chains.RData")
load("arthritis_10chains.RData")
source('R/aaa_generics.R') 
source('R/compute_LL.R') 
source('R/func_aux.R') 
source('R/func_distcond.R') 
source('R/func_FIM.R')
source('R/func_plots.R') 
source('R/func_simulations.R') 
source('R/main.R')
source('R/main_estep.R')
source('R/main_estep_incremental.R')
source('R/main_initialiseMainAlgo.R') 
source('R/main_mstep.R') 
source('R/SaemixData.R')
source('R/SaemixModel.R') 
source('R/SaemixRes.R') 
source('R/SaemixObject.R') 
source('R/zzz.R') 

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

project.file <- "mlxProjects/arthritis/arthritis_projet.mlxtran"
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

arthritis_data <- readDatamlx(project = project.file)
treat <- arthritis_data$treatment[,c(3)]
# arthritis.saemix <- cbind(arthritis_data$Y,treat)
arthritis.saemix <- arthritis_data$y

saemix.data<-saemixData(name.data=arthritis.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("time"),name.response=c("y"), name.X="time")

saemix.model<-saemixModel(model=model1cpt,description="arthritis",type="structural"
  ,psi0=matrix(c(10,10),ncol=2,byrow=TRUE, dimnames=list(NULL, c("tauD","tauS"))),
  fixed.estim=c(1,1),transform.par=c(1,1),omega.init=matrix(diag(2),ncol=2,byrow=TRUE),
  covariance.model=matrix(diag(2),ncol=2,byrow=TRUE), error.model = "proportional")

K1 = 1000
K2 = 500
iterations = 1:(K1+K2)
end = K1+K2

runtime = 20
nchains = 1
options<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=100,sampling='randompass', duration = runtime)
arthritis<-data.frame(saemix(saemix.model,saemix.data,options))
arthritis <- cbind(iterations, arthritis[-1,])
row_sub_ref  = apply(arthritis, 1, function(row) all(row !=0 ))
arthritis <- arthritis[row_sub_ref,]
arthritis$algo <- 'full'
arthritis$iterations <- seq(0,runtime, length.out=length(arthritis$iterations))


nchains = 1
options50<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=50,sampling='randompass', duration = runtime)
arthritis50<-data.frame(saemix(saemix.model,saemix.data,options50))
arthritis50 <- cbind(iterations, arthritis50[-1,])
row_sub_50  = apply(arthritis50, 1, function(row) all(row !=0 ))
arthritis50 <- arthritis50[row_sub_50,]
arthritis50$algo <- 'half'
arthritis50$iterations <- seq(0,runtime, length.out=length(arthritis50$iterations))


options25<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
 nb.replacement=25,sampling='randompass', duration = runtime)
arthritis25<-data.frame(saemix(saemix.model,saemix.data,options25))
arthritis25 <- cbind(iterations, arthritis25[-1,])
row_sub_25  = apply(arthritis25, 1, function(row) all(row !=0 ))
arthritis25 <- arthritis25[row_sub_25,]
arthritis25$algo <- 'quarter'
arthritis25$iterations <- seq(0,runtime, length.out=length(arthritis25$iterations))




# +scale_x_log10()
dim(arthritis)
dim(arthritis50)
dim(arthritis25)

comparison <- 0
# comparison <- rbind(arthritis,arthritis)
# comparison <- rbind(arthritis,arthritis50)
load("arthritis_10chains.RData")
load("arthritis_40chains.RData")
seplot <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  graf <- ggplot(df)+geom_line(aes(iterations,value,by=value,colour = df$algo),show.legend = legend) +
  xlab("iterations") +scale_x_log10()+ ylab('value') + facet_wrap(~variable,scales = "free_y") + theme_bw()
  grid.arrange(graf)
  # do.call("grid.arrange", c(graf, ncol=1, top=title))
}
comparison <- rbind(arthritis,arthritis25,arthritis50)
var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
prec <- seplot(var, title="comparison",legend=TRUE)


write.csv(comparison, file = "../notebooks/arthritis.csv")
