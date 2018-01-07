setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/averageiterate")
source('mixtureFunctions.R') 

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/averageiterate/avg")
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
  source('main_initialiseMainAlgoavg.R')
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 

library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
###RTTE
timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/paramToRV/data/rtte_data.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
saemix.data_rtte<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
timetoevent.model<-function(psi,id,xidep) {
  T<-xidep[,1]
  y<-xidep[,2]
  N <- nrow(psi)
  Nj <- length(T)
  censoringtime = 20
  lambda <- psi[id,1]
  beta <- psi[id,2]
  init <- which(T==0)
  cens <- which(T==censoringtime)
  ind <- setdiff(1:Nj, append(init,cens))
  hazard <- (beta/lambda)*(T/lambda)^(beta-1)
  H <- (T/lambda)^beta
  logpdf <- rep(0,Nj)
  logpdf[cens] <- -H[cens] + H[cens-1]
  logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])
  return(logpdf)
}

saemix.model_rtte<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))

##RUNS

K1 = 100
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2


#With var no sa
options.ref<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=0)
rtte.ref<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options.ref))
rtte.ref <- cbind(iterations, rtte.ref)
graphConvMC_twokernels(rtte.ref,rtte.ref)

options.avg<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=1)
rtte.avg<-data.frame(saemix(saemix.model_rtte,saemix.data_rtte,options.avg))
rtte.avg <- cbind(iterations, rtte.avg)

graphConvMC_twokernels(rtte.avg,rtte.avg)

final_rwm <- 0
var_rwm <- 0
error_rwm <- 0
lambda_true = 10
o_lambda_true = 1
beta_true = 2
o_beta_true = 1
final_mix <- 0
true_param <- c(lambda_true,beta_true,o_lambda_true,o_beta_true)
var_mix <- 0
error_mix <- 0


seed0 = 39546
replicate = 4

for (j in 1:replicate){

  model2 <- inlineModel("

  [LONGITUDINAL]
  input = {beta,lambda}  

  EQUATION:
  h=(beta/lambda)*(t/lambda)^(beta-1)

  DEFINITION:
  e = {type               = event, 
       rightCensoringTime = 20,  
       hazard             = h}
  [INDIVIDUAL]
  input={lambda_pop, o_lambda,beta_pop, o_beta}
                        
  DEFINITION:
  lambda  ={distribution=lognormal, prediction=lambda_pop,  sd=o_lambda}
  beta  ={distribution=lognormal, prediction=beta_pop,  sd=o_beta}
       ")


  p <- c(lambda_pop=4, o_lambda= 0.2,
         beta_pop = 2,o_beta = 0.3)
  h <- list(name='h', time=seq(0, 6, by=1))
  e <- list(name='e', time=0)

  N <- 200
  res <- simulx(model     = model2, 
                settings  = list(seed=j*123),
                parameter = p, 
                output    = list(h,e), 
                 group     = list(size = N))

res
writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Desktop/CSDA_code/rtte/rtte2.csv")
head(read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code/rtte/rtte2.csv", header=T, sep=","))
timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code/rtte/rtte2.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
  saemix.data<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
  print(j)
  options.ref<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,map.range=c(0),nbiter.burn =0, av=0,avg=0)
  rtte.ref<-data.frame(saemix(saemix.model_rtte,saemix.data,options.ref))
  rtte.ref <- cbind(iterations, rtte.ref)
  rtte.ref['individual'] <- j
  final_rwm <- rbind(final_rwm,rtte.ref)

  var_rwm <- var_rwm + (rtte.ref[,2:5]-true_param)^2
  ML <- rtte.ref[,2:5]
  ML[1:(end+1),]<- rtte.ref[end+1,2:5]
  error_rwm <- error_rwm + (rtte.ref[,2:5]-ML)^2
  

  # options.avg<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(0), av=0,avg=1)
  # rtte.avg<-data.frame(saemix(saemix.model_rtte,saemix.data,options.avg))
  # rtte.avg <- cbind(iterations, rtte.avg)
  # rtte.avg['individual'] <- j
  # final_mix <- rbind(final_mix,rtte.avg)
  # var_mix <- var_mix + (rtte.avg[,2:5]-true_param)^2
  # ML <- rtte.avg[,2:5]
  # ML[1:(end+1),]<- rtte.avg[end+1,2:5]
  # error_mix <- error_mix + (rtte.avg[,2:5]-ML)^2
  
}

graphConvMC_twokernels(rtte.ref,rtte.ref)

first <- subset(final_rwm, individual==1)
second <- subset(final_rwm, individual==2)

first[150:end,] - second[150:end,]


graphConvMC_diff(final_rwm,final_rwm)








error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix

error_rwm <- cbind(iterations, error_rwm)
error_mix <- cbind(iterations, error_mix)

err_mix<- rtte.ref
err_rwm<- rtte.ref
err_rwm[,2:5] <- error_rwm[,2:5]
err_mix[,2:5] <- error_mix[,2:5]

err_mix[2,] = err_rwm[2,]

graphConvMC_se1 <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=2) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=2)+
      xlab("") + ylab(expression(paste(lambda)))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=15, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=30)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}

graphConvMC_se2 <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=2) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=2)+
      xlab("") + ylab(expression(paste(omega,".",lambda)))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=15, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=30)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}

# graphConvMC_diff2(err_rwm[-1,c(1,2,4,6)],err_mix[-1,c(1,2,4,6)], title="Quadratic errors RTTE")
# graphConvMC_diff2(err_rwm[-1,c(1,2,4,6)],err_mix[-1,c(1,2,4,6)])


c <- graphConvMC_se1(err_rwm[K1:end,c(1,2,6)],err_mix[K1:end,c(1,2,6)])
d <- graphConvMC_se2(err_rwm[K1:end,c(1,4,6)],err_mix[K1:end,c(1,4,6)])

grid.arrange(c,d, ncol=2)


c <- graphConvMC_se1(err_rwm[-1,c(1,2,6)],err_mix[-1,c(1,2,6)])
d <- graphConvMC_se2(err_rwm[-1,c(1,4,6)],err_mix[-1,c(1,4,6)])

grid.arrange(c,d, ncol=2)





