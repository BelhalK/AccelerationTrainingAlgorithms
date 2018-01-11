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
  source('main_avg.R')
  source('main_estep.R')
  source('main_initialiseMainAlgo.R')
  source('main_initialiseMainAlgoavg.R')
  source('main_mstep.R') 
  source('main_mstep_avg.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 


library("mlxR")
library(sgd)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
###WARFA
warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/averageiterate/data/warfarin_data.txt", header=T)
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
  ,psi0=matrix(c(1,4,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))




K1 = 100
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2


#With var no sa
options.ref<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=0)
warfa.ref<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options.ref))
warfa.ref <- cbind(iterations, warfa.ref)
graphConvMC_twokernels(warfa.ref,warfa.ref)

options.avg<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=1)
warfa.avg<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options.avg))
warfa.avg <- cbind(iterations, warfa.avg)


var_rwm <- 0
error_rwm <- 0


var_mix <- 0
error_mix <- 0

var_avgsa <- 0
error_avgsa <- 0

final_rwm <- 0
final_mix <- 0
final_avgsa <- 0

replicate <- 30
for (m in 1:replicate){
  setwd("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin")
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

  setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/averageiterate/avg")
  warfarin.saemix_less <- warfarin.saemix[,]
  saemix.data_warfa<-saemixData(name.data=warfarin.saemix_less,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")


  options.ref<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1,nbiter.mcmc = c(2,2,2,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=0, map.range=c(0))
  warfa.ref<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options.ref)$parpop)
  warfa.ref <- cbind(iterations, warfa.ref)
  # var_rwm <- var_rwm + (warfa.ref[,2:8]-true_param)^2
  ML <- warfa.ref[,2:8]
  ML[1:(end+1),]<- warfa.ref[end+1,2:8]
  error_rwm <- error_rwm + (warfa.ref[,2:8]-ML)^2
  warfa.ref['individual'] <- j
  final_rwm <- rbind(final_rwm,warfa.ref)
  
  options.avgsa<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1,nbiter.mcmc = c(2,2,2,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=1, map.range=c(0))
  warfa.avgsa<-data.frame(saemix_avg(saemix.model_warfa,saemix.data_warfa,options.avgsa))
  warfa.avgsa <- cbind(iterations, warfa.avgsa)
  # var_mix <- var_mix + (theo_mix[,2:8]-true_param)^2
  ML <- warfa.avgsa[,2:8]
  ML[1:(end+1),]<- warfa.avgsa[end+1,2:8]
  error_avgsa <- error_avgsa + (warfa.avgsa[,2:8]-ML)^2
  warfa.avgsa['individual'] <- m
  final_avgsa <- rbind(final_avgsa,warfa.avgsa)

  options.avg<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1,nbiter.mcmc = c(2,2,2,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, av=0,avg=1, map.range=c(0))
  warfa.avg<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options.avg)$newparpop)
  warfa.avg <- cbind(iterations, warfa.avg)
  # var_mix <- var_mix + (theo_mix[,2:8]-true_param)^2
  ML <- warfa.avg[,2:8]
  ML[1:(end+1),]<- warfa.avg[end+1,2:8]
  error_mix <- error_mix + (warfa.avg[,2:8]-ML)^2
  warfa.avg['individual'] <- m
  final_mix <- rbind(final_mix,warfa.avg)

  print(warfa.ref[50:K1,2:8] - warfa.avg[50:K1,2:8])

}
ML_rwm <- subset(final_rwm, iterations=200)
ML_avg <- subset(final_mix, iterations=200)
ML_avgsa <- subset(final_avgsa, iterations=200)

graphConvMC_twokernels(warfa.ref,warfa.avg)

error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix
error_avgsa <- 1/replicate*error_avgsa

error_rwm <- cbind(iterations, error_rwm)
error_mix <- cbind(iterations, error_mix)
error_avgsa <- cbind(iterations, error_avgsa)

err_mix<- warfa.ref
err_rwm<- warfa.ref
err_avgsa<- warfa.ref

err_rwm[,2:8] <- error_rwm[,2:8]
err_mix[,2:8] <- error_mix[,2:8]
err_avgsa[,2:8] <- error_avgsa[,2:8]

err_mix[2,] = err_rwm[2,] = err_avgsa[2,]


a <- err_rwm[-1,]
b <- err_mix[-1,]
avgsa <- err_avgsa[-1,]


graphConvMC_diff3 <- function(df,df2, title=NULL, ylim=NULL)
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
      xlab("") + ylab("")  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=14, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=14, angle=0))+ggtitle(names(df[j]))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


graphConvMC_diff4 <- function(df,df2, title=NULL, ylim=NULL)
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
      xlab("") +scale_y_continuous( limits=c(2, 8))+ ylab("")  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=14, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=14, angle=0))+ggtitle(names(df[j]))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}

graphConvMC_diff5 <- function(df,df2,df3, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  df3$individual <- as.factor(df3$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=2) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=2)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="green",linetype = 2,size=2)+
      xlab("") + ylab("")  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=14, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=14, angle=0))+ggtitle(names(df[j]))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}

c <- graphConvMC_diff3(a[K1:end,c(1,3,9)],b[K1:end,c(1,3,9)])
d <- graphConvMC_diff3(a[K1:end,c(1,6,9)],b[K1:end,c(1,6,9)])
e <- graphConvMC_diff5(a[K1:end,c(1,3,9)],b[K1:end,c(1,3,9)],avgsa[K1:end,c(1,3,9)])
f <- graphConvMC_diff5(a[K1:end,c(1,6,9)],b[K1:end,c(1,6,9)],avgsa[K1:end,c(1,6,9)])

# grid.arrange(c,d, ncol=2)
grid.arrange(e,f, ncol=2)

c <- graphConvMC_diff3(a[(K1-30):end,c(1,3,9)],b[(K1-30):end,c(1,3,9)])
d <- graphConvMC_diff3(a[(K1-30):end,c(1,6,9)],b[(K1-30):end,c(1,6,9)])
e <- graphConvMC_diff5(a[(K1-30):end,c(1,3,9)],b[(K1-30):end,c(1,3,9)],avgsa[(K1-30):end,c(1,3,9)])
f <- graphConvMC_diff5(a[(K1-30):end,c(1,6,9)],b[(K1-30):end,c(1,6,9)],avgsa[(K1-30):end,c(1,6,9)])
grid.arrange(c,d, ncol=2)


c <- graphConvMC_diff3(a[1:end,c(1,3,9)],b[1:end,c(1,3,9)])
d <- graphConvMC_diff3(a[1:end,c(1,6,9)],b[1:end,c(1,6,9)])
e <- graphConvMC_diff5(a[1:end,c(1,3,9)],b[1:end,c(1,3,9)],avgsa[1:end,c(1,3,9)])
f <- graphConvMC_diff5(a[1:end,c(1,6,9)],b[1:end,c(1,6,9)],avgsa[1:end,c(1,6,9)])
grid.arrange(c,d, ncol=2)