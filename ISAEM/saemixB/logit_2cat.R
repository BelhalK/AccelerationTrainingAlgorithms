
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
  source('mixtureFunctions.R')

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/")
source("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/plots.R")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################

catModel <- inlineModel("
[LONGITUDINAL]
input =  {beta0,gamma0,delta0, dose}
dose = {use=regressor}
EQUATION:
lm0 = beta0+gamma0*t + delta0*dose

D = exp(lm0)+1

p0 = exp(lm0)/D
p1 = 1/D

DEFINITION:
y = {type=categorical, categories={0, 1}, 
     P(y=0)=p0,
     P(y=1)=p1}

[INDIVIDUAL]
input={beta0_pop, o_beta0,
      gamma0_pop, o_gamma0,
      delta0_pop, o_delta0}


DEFINITION:
beta0  ={distribution=normal, prediction=beta0_pop,  sd=o_beta0}

gamma0  ={distribution=normal, prediction=gamma0_pop,  sd=o_gamma0}

delta0  ={distribution=normal, prediction=delta0_pop,  sd=o_delta0}

")


nobs = 10
t<- seq(1, 100, by=nobs)
reg <- list(name='dose',
            time=t,
            value=20)

out  <- list(name='y', time=seq(1, 100, by=nobs))
N  <- 1000
p <- c(beta0_pop=1, o_beta0=0, 
       gamma0_pop= -1, o_gamma0=0.5,
       delta0_pop=3, o_delta0=0.4)

g1 <- list(size=N, parameter=p)
res <- simulx(model=catModel, regressor = reg, output=out, group=g1)
plot1 <- catplotmlx(res$y)
print(plot1)


writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/cat_data.csv")
cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/cat_data.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"), name.predictors=c("y","dose","time"))



cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]
dose<-xidep[,2]
time<-xidep[,3]

beta0 <- psi[id,1]

gamma0 <- psi[id,2]

delta0 <- psi[id,3]

lm0 <- beta0+gamma0*time + delta0*dose

D <- exp(lm0)+1

P0 <- exp(lm0)/D
P1 <- 1/D

P.obs = (level==0)*P0+(level==1)*P1

return(P.obs)
}


cov <- matrix(c(1,0,0,
                0,1,0,
                0,0,1),ncol=3, byrow=TRUE)

saemix.model<-saemixModel(model=cat_data.model,description="cat model",type="likelihood",   
  psi0=matrix(c(1,-2,1),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("beta0",
    "gamma0",
    "delta0"))), 
  transform.par=c(0,0,0), fixed.estim=c(1,1,1),covariance.model=cov,omega.init=cov,error.model="constant")



K1 = 200
K2 = 1

iterations = 1:(K1+K2+1)
end = K1+K2
seed0 = 444


options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=100,sampling='randompass')
theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

options.incremental<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
theo_mix <- cbind(iterations, theo_mix)

options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-data.frame(saemix(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_ref)
 

graphConvMC_twokernels(theo_ref,theo_mix)

replicate = 1

final_rwm <- 0
final_incremental <- 0
final_incremental25 <- 0
for (m in 1:replicate){
  print(m)
  print(m)
  # l = list(c(1,5,1,0,0,0),c(0.8,12,0.8,0,0,0),c(1.2,3,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  l = list(c(1,-2,1),c(3,3,3),c(1,2.5,2.5))
  
saemix.model<-saemixModel(model=cat_data.model,description="cat model",   
  psi0=matrix(c(1,-2,1),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("beta0",
    "gamma0",
    "delta0"))), 
  transform.par=c(0,0,0), fixed.estim=c(0,1,1),covariance.model=cov,omega.init=cov,error.model="constant")




  options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=100,sampling='randompass')
  theo_ref<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- m
  theo_ref <- theo_ref[-1,]
  theo_ref_scaled <- theo_ref
  theo_ref_scaled$iterations = seq(1, 4*end, by=4)
  theo_ref_scaled <- theo_ref_scaled[rep(seq_len(nrow(theo_ref_scaled)), each=4),]

  final_rwm <- rbind(final_rwm,theo_ref_scaled[0:end,])

  options.incremental<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=50,sampling='randompass')
  theo_mix<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental))
  theo_mix <- cbind(iterations, theo_mix)
  theo_mix['individual'] <- m
  theo_mix <- theo_mix[-1,]
  theo_mix_scaled <- theo_mix
  theo_mix_scaled$iterations = seq(1, 2*end, by=2)
  theo_mix_scaled <- theo_mix_scaled[rep(seq_len(nrow(theo_mix_scaled)), each=2),]
  final_incremental <- rbind(final_incremental,theo_mix_scaled[0:end,])

  options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
  theo_mix25<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental25))
  theo_mix25 <- cbind(iterations, theo_mix25)
  theo_mix25['individual'] <- m
  theo_mix25 <- theo_mix25[-1,]
  theo_mix25$iterations = seq(1, end, by=1)
  theo_mix25[end,] <- theo_mix25[(end-2),]
  final_incremental25 <- rbind(final_incremental25,theo_mix25[0:end,])
}


a <- graphConvMC_diffz(final_rwm[,c(1,4,7)],final_incremental[,c(1,4,7)],final_incremental25[,c(1,4,7)])

#MC STUDY

graphConvMC_sec <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype=1,size=1)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="green",linetype=1,size=1)+
      xlab("") + scale_x_log10()+ylab(expression(paste(beta,"2")))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                          size=30, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                          size=30, angle=0))+theme(axis.title = element_text(color="black", face="bold", size=30)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


graphConvMC_sed <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype=1,size=1)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="green",linetype=1,size=1)+
      xlab("") + scale_x_log10()+ ylab(expression(paste(omega,"2")))   + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                          size=30, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                          size=30, angle=0))+theme(axis.title = element_text(color="black", face="bold", size=30)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}

K1 = 600
K2 = 0

iterations = 1:(K1+K2+1)
end = K1+K2

final_rwm <- 0
final_ref <- 0
error_rwm <- 0

th1 <- 1
th2 <- 3
th3 <- 1
o_th1 <- 1
o_th2 <- 0.5
o_th3 <- 0.5
final_mix <- 0
final_mix25 <- 0
final_mix10 <- 0
true_param <- c(th1,th2,th3,o_th1,o_th2,o_th3)
error_mix <- 0
error_mix25 <- 0
error_mix10 <- 0


seed0 = 39546
replicate = 20

for (j in 1:replicate){

catModel <- inlineModel("
[LONGITUDINAL]
input =  {beta0,gamma0,delta0, dose}
dose = {use=regressor}
EQUATION:
lm0 = beta0+gamma0*t + delta0*dose

D = exp(lm0)+1

p0 = exp(lm0)/D
p1 = 1/D

DEFINITION:
y = {type=categorical, categories={0, 1}, 
     P(y=0)=p0,
     P(y=1)=p1}

[INDIVIDUAL]
input={beta0_pop, o_beta0,
      gamma0_pop, o_gamma0,
      delta0_pop, o_delta0}


DEFINITION:
beta0  ={distribution=normal, prediction=beta0_pop,  sd=o_beta0}
gamma0  ={distribution=normal, prediction=gamma0_pop,  sd=o_gamma0}
delta0  ={distribution=normal, prediction=delta0_pop,  sd=o_delta0}

")


nobs = 10
t<- seq(1, 100, by=nobs)
reg <- list(name='dose',
            time=t,
            value=20)

out  <- list(name='y', time=seq(1, 100, by=nobs))
N  <- 10000
p <- c(beta0_pop=1, o_beta0=0.3, 
       gamma0_pop= -1, o_gamma0=0.5,
       delta0_pop=3, o_delta0=0.4)

g1 <- list(size=N, parameter=p)
res <- simulx(model=catModel, regressor = reg, output=out, group=g1)


writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/logistic_2cat_se.csv")
cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/logistic_2cat_se.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"), name.predictors=c("y","dose","time"))



cov <- matrix(c(0,0,0,
                0,1,0,
                0,0,1),ncol=3, byrow=TRUE)

saemix.model<-saemixModel(model=cat_data.model,description="cat model",   
  psi0=matrix(c(2,-2,1),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("beta0",
    "gamma0",
    "delta0"))), 
  transform.par=c(0,0,0), fixed.estim=c(0,1,1),covariance.model=cov,omega.init=cov,error.model="constant")

  print(j)
  options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=100,sampling='randompass')
  theo_ref<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  ML <- theo_ref[,2:6]
  ML[1:(end+1),]<- theo_ref[end,2:6]
  error_rwm <- error_rwm + (theo_ref[,2:6]-ML)^2
  theo_ref['individual'] <- j
  final_ref <- rbind(final_ref,theo_ref)


  options.incremental<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=50,sampling='randompass')
  theo_mix<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental))
  theo_mix <- cbind(iterations, theo_mix)
  ML <- theo_mix[,2:6]
  ML[1:(end+1),]<- theo_mix[end,2:6]
  error_mix <- error_mix + (theo_mix[,2:6]-ML)^2
  theo_mix['individual'] <- j
  final_mix <- rbind(final_mix,theo_mix)
  
  options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
  theo_mix25<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental25))
  theo_mix25 <- cbind(iterations, theo_mix25)
  ML <- theo_mix25[,2:6]
  ML[1:(end+1),]<- theo_mix25[end,2:6]
  error_mix25 <- error_mix25 + (theo_mix25[,2:6]-ML)^2
  theo_mix25['individual'] <- j
  final_mix25 <- rbind(final_mix25,theo_mix25)

  # options.incremental10<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=10)
  # theo_mix10<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental10))
  # theo_mix10 <- cbind(iterations, theo_mix10)
  # ML <- theo_mix10[,2:6]
  # ML[1:(end+1),]<- theo_mix10[end+1,2:6]
  # error_mix10 <- error_mix10 + (theo_mix10[,2:6]-ML)^2
  # theo_mix10['individual'] <- j
  # final_mix10 <- rbind(final_mix10,theo_mix10)
}

error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix
error_mix25 <- 1/replicate*error_mix25
# error_mix10 <- 1/replicate*error_mix10


err_mix<- theo_ref[-1,]
err_rwm<- theo_ref[-1,]
err_mix25<- theo_ref[-1,]

err_rwm[,2:6] <- error_rwm[-1,1:5]
err_mix[,2:6] <- error_mix[-1,1:5]
err_mix25[,2:6] <- error_mix25[-1,1:5]


err_rwm_scaled <- err_rwm
err_rwm_scaled$iterations = seq(1, 4*end, by=4)
err_mix_scaled <- err_mix
err_mix_scaled$iterations = seq(1, 2*end, by=2)
err_mix25$iterations = 1:((K1+K2))



err_mix_scaled[1,] = err_mix25[1,] = err_rwm_scaled[1,]



for (i in 2:6){
# i = 6
prec <- graphConvMC_sec_icml(err_rwm_scaled[0:(end-1),c(1,i,7)],err_mix_scaled[0:(end-1),c(1,i,7)],err_mix25[0:(end-1),c(1,i,7)])
# assign(paste("prec", i, sep = ""), prec) 
setwd("/Users/karimimohammedbelhal/Desktop/")
ggsave(paste("logit2cat_", i, ".png", sep=""),prec)
}


prec1 <- graphConvMC_sec_icml(err_rwm_scaled[0:(end-1),c(1,3,7)],err_mix_scaled[0:(end-1),c(1,3,7)],err_mix25[0:(end-1),c(1,3,7)])
prec2 <- graphConvMC_sec_icml(err_rwm_scaled[0:(end-1),c(1,4,7)],err_mix_scaled[0:(end-1),c(1,4,7)],err_mix25[0:(end-1),c(1,4,7)])
prec3 <- graphConvMC_sec_icml(err_rwm_scaled[0:(end-1),c(1,5,7)],err_mix_scaled[0:(end-1),c(1,5,7)],err_mix25[0:(end-1),c(1,5,7)])
prec4 <- graphConvMC_sec_icml(err_rwm_scaled[0:(end-1),c(1,6,7)],err_mix_scaled[0:(end-1),c(1,6,7)],err_mix25[0:(end-1),c(1,6,7)])

quad <- grid.arrange(prec1,prec2,prec3,prec4,ncol=2)
ggsave(paste("logit2cat_quadrant.png", sep=""),quad)

graphConvMC_sec_icml <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",linetype= "solid",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="black",linetype="longdash",size=1)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype="dotted",size=1)+
      xlab("") + scale_x_log10()+ylab(names(df[j]))   + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=15, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=15)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


graphConvMC_sed_icml <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",linetype= "solid",size=2) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="black",linetype="longdash",size=2)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype="dotted",size=2)+
      xlab("") + scale_x_log10()+ ylab(expression(paste(omega,"2")))   + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=30, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=30, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=30)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
    do.call("grid.arrange", c(graf, ncol=1, top=title))
}
