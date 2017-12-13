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
source('main_time.R')
  source('main_estep_time.R')
  source('main_mstep_time.R') 
  source('func_aux_time.R') 
  source('SaemixObject_time.R') 
  source('main_initialiseMainAlgo_time.R') 
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



library(saemix)
PD1.saemix<-read.table( "data/PD1.saemix.tab",header=T,na=".")
saemix.data1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))

PD2.saemix<-read.table( "data/PD2.saemix.tab",header=T,na=".")
PD2.saemix <- PD2.saemix[1:168,]
saemix.data2<-saemixData(name.data=PD2.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))


PD3.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/pd1.csv", header=T, sep=",")
saemix.data3<-saemixData(name.data=PD3.saemix,header=TRUE,name.group=c("id"),
name.predictors=c("time"),name.response=c("y"),
units=list(x="mg",y="-",covariates="-"))

modelemax<-function(psi,id,xidep) {
# input:
# psi : matrix of parameters (3 columns, E0, Emax, EC50)
# id : vector of indices
# xidep : dependent variables (same nb of rows as length of id)
# returns:
# a vector of predictions of length equal to length of id
dose<-xidep[,1]
e0<-psi[id,1]
emax<-psi[id,2]
e50<-psi[id,3]
f<-e0+emax*dose/(e50+dose)
return(f)
}

saemix.model<-saemixModel(model=modelemax,description="Emax model",
psi0=matrix(c(20,300,25,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")


saemix.model1<-saemixModel(model=modelemax,description="Emax model",
psi0=matrix(c(10,200,25,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")

K1 = 100
K2 = 50
iterations = 1:(K1+K2+1)
gd_step = 0.01

end = K1+K2

seed0 = 39546
#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,0,0,0,0,0),displayProgress=FALSE, nbiter.saemix = c(K1,K2),nbiter.burn =0)
theo_ref<-data.frame(saemix_mamyula(saemix.model,saemix.data1,options))
theo_ref <- cbind(iterations, theo_ref)


options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,0,0,0,0,0),displayProgress=FALSE, nbiter.saemix = c(K1,K2),nbiter.burn =0)
theo_ref<-data.frame(saemix_mamyula(saemix.model1,saemix.data3,options))
theo_ref <- cbind(iterations, theo_ref)

theo_ref[end,]

graphConvMC_twokernels(theo_ref,theo_ref, title="new kernel")
#saem with mala
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,2,0,0),displayProgress=FALSE,nbiter.saemix = c(K1,K2),sigma.val = 0.01,gamma.val=0.01,nbiter.burn =0)
theo_mala<-data.frame(saemix_mamyula(saemix.model,saemix.data1,options.mala))
theo_mala <- cbind(iterations, theo_mala)


#saem with mamyula
options.mamyula<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,0,2,0),displayProgress=FALSE,nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01,lambda.val=0.2,nbiter.burn =0)
theo_mamyula<-data.frame(saemix_mamyula(saemix.model,saemix.data1,options.mamyula))
theo_mamyula <- cbind(iterations, theo_mamyula)

options.mamyula<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,0,2,0),displayProgress=FALSE,nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01,lambda.val=0.2,nbiter.burn =0)
theo_mamyula<-data.frame(saemix_mamyula(saemix.model1,saemix.data3,options.mamyula))
theo_mamyula <- cbind(iterations, theo_mamyula)


theo_ref[end,]
theo_mamyula[end,]

graphConvMC_twokernels(theo_ref,theo_mala, title="new kernel")
graphConvMC_threekernels(theo_ref,theo_mala,theo_mamyula, title="new kernel")
graphConvMC_threekernels(theo_ref,theo_mamyula,theo_mamyula, title="new kernel")



replicate = 3

final_rwm <- 0
final_mix <- 0
for (m in 1:replicate){
  print(m)
  print(m)
  l = list(c(15,250,20,0,0,0),c(10,200,25,0,0,0),c(5,170,15,0,0,0),c(20,300,20,0,0,0))
  

saemix.model1<-saemixModel(model=modelemax,description="Emax model",
psi0=matrix(l[[m]],ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")


  options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,0,0,0,0,0),displayProgress=FALSE, nbiter.saemix = c(K1,K2),nbiter.burn =0)
theo_ref<-data.frame(saemix_mamyula(saemix.model1,saemix.data3,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,theo_ref[-1,])

  options.mamyula<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,0,2,0),displayProgress=FALSE,nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01,lambda.val=0.2,nbiter.burn =0)
theo_mamyula<-data.frame(saemix_mamyula(saemix.model1,saemix.data3,options.mamyula))
  theo_mix <- cbind(iterations, theo_mamyula)
  theo_mix['individual'] <- m
  final_mix <- rbind(final_mix,theo_mix[-1,])
}

graphConvMC_diff2(final_rwm[,],final_mix[,], title="Diff intial param pd")





a <- graphConvMC_diff4(final_rwm[,c(1,2,9)],final_mix[,c(1,2,9)])
b <- graphConvMC_diff3(final_rwm[,c(1,5,9)],final_mix[,c(1,5,9)])

grid.arrange(a,b, ncol=2)



replicate = 3

final_rwm <- 0
final_mix <- 0
for (m in 1:replicate){
  print(m)
  print(m)
  l = list(c(15,280,15,0,0,0),c(25,300,25,0,0,0),c(20,320,20,0,0,0),c(20,300,20,0,0,0))
  saemix.model<-saemixModel(model=modelemax,description="Emax model",
psi0=matrix(l[[m]],ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")


  options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,0,0,0,0,0),displayProgress=FALSE, nbiter.saemix = c(K1,K2),nbiter.burn =0)
theo_ref<-data.frame(saemix_mamyula(saemix.model,saemix.data1,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,theo_ref[-1,])

  options.mamyula<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,0,2,0),displayProgress=FALSE,nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01,lambda.val=0.2,nbiter.burn =0)
theo_mamyula<-data.frame(saemix_mamyula(saemix.model,saemix.data1,options.mamyula))
  theo_mix <- cbind(iterations, theo_mamyula)
  theo_mix['individual'] <- m
  final_mix <- rbind(final_mix,theo_mix[-1,])
}



graphConvMC_diff2(final_rwm[,c(1,4,9)],final_mix[,c(1,4,9)], title="Diff intial param pd")




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
      xlab("") +scale_x_log10()+ ylab("w2.V")  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=14, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=14, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=22)) 
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
      xlab("") +scale_x_log10()+scale_y_continuous()+ ylab(names(df[j]))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=14, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=14, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=22)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


a <- graphConvMC_diff4(final_rwm[,c(1,4,9)],final_mix[,c(1,4,9)])
b <- graphConvMC_diff3(final_rwm[,c(1,7,9)],final_mix[,c(1,7,9)])

grid.arrange(a,b, ncol=2)

#run on diff datasets to compute residual errors

var_rwm <- 0
error_rwm <- 0


var_mix <- 0
error_mix <- 0

replicate = 40

final_rwm <- 0
final_mix <- 0

replicate <- 40
for (m in 1:replicate){
  


model2 <- inlineModel("
                      [LONGITUDINAL]
                      input = {e, em, ec,a}

                      EQUATION:
                      Cc = e+(em*t)/(ec+t)
                      
                      DEFINITION:
                      y1 ={distribution=normal, prediction=Cc, sd=a}
                      
                      [INDIVIDUAL]
                      input={e_pop,o_e,em_pop,o_em,ec_pop,o_ec}
                      
                      DEFINITION:
                      e   ={distribution=lognormal, prediction=e_pop,   sd=o_e}
                      em   ={distribution=lognormal, prediction=em_pop,   sd=o_em}
                      ec  ={distribution=lognormal, prediction=ec_pop,  sd=o_ec}
                      ")

adm  <- list(amount=c(0, 10, 90), time=c(0, 10, 90))


p <- c(e_pop=20, o_e=0.5,
       em_pop=100, o_em=0.5, 
       ec_pop=10, o_ec=0.5,  
       a=0.1)
y1 <- list(name='y1', time=c(0, 10, 90))


res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=50, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/pdse.csv")


#modification for mlxsaem dataread function
obj <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/pdse.csv", header=T, sep=",")
a <- which(obj[,4]==0)
data <- obj[-a,]
b <- which(data[,4]==10)
data <- data[-b,]
c <- which(data[,4]==90)
data <- data[-c,]

write.table(data, "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/pdse.csv", sep=",", row.names=FALSE,quote = FALSE, col.names=TRUE)


PD4.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/pdse.csv", header=T, sep=",")
saemix.data4<-saemixData(name.data=PD4.saemix,header=TRUE,name.group=c("id"),
name.predictors=c("time"),name.response=c("y"),
units=list(x="mg",y="-",covariates="-"))

saemix.modelse<-saemixModel(model=modelemax,description="Emax model",
psi0=matrix(c(10,200,25,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")





  options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,0,0,0,0,0),displayProgress=FALSE, nbiter.saemix = c(K1,K2),nbiter.burn =0)
  theo_ref<-data.frame(saemix_mamyula(saemix.modelse,saemix.data4,options))
  theo_ref <- cbind(iterations, theo_ref)
  # var_rwm <- var_rwm + (theo_ref[,2:8]-true_param)^2
  ML <- theo_ref[,2:8]
  ML[1:(end+1),]<- theo_ref[end+1,2:8]
  error_rwm <- error_rwm + (theo_ref[,2:8]-ML)^2
  theo_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,theo_ref)
  
  
  options.mamyula<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,0,2,0),displayProgress=FALSE,nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01,lambda.val=0.2,nbiter.burn =0)
  print(m)
  theo_mix<-data.frame(saemix_mamyula(saemix.modelse,saemix.data4,options.mamyula))
  theo_mix <- cbind(iterations, theo_mix)
  # var_mix <- var_mix + (theo_mix[,2:8]-true_param)^2
  ML <- theo_mix[,2:8]
  ML[1:(end+1),]<- theo_mix[end+1,2:8]
  error_mix <- error_mix + (theo_mix[,2:8]-ML)^2
  theo_mix['individual'] <- m
  final_mix <- rbind(final_mix,theo_mix)
}

error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix

error_rwm <- cbind(iterations, error_rwm)
error_mix <- cbind(iterations, error_mix)

err_mix<- theo_ref
err_rwm<- theo_ref
err_rwm[,2:8] <- error_rwm[,2:8]
err_mix[,2:8] <- error_mix[,2:8]

err_mix[2,] = err_rwm[2,]


# graphConvMC_diff(err_rwm[-1,],err_mix[-1,], title="Quadratic errors Warfa")
# graphConvMC_diff2(err_rwm[-1,],err_mix[-1,], title="Quadratic errors Warfa")

# graphConvMC_diff2(err_rwm[-1,c(1,3,6,9)],err_mix[-1,c(1,3,6,9)])

a <- err_rwm[-1,]
b <- err_mix[-1,]
# graphConvMC_diff2(err_rwm[-1,c(1,3,6,9)],err_mix[-1,c(1,3,6,9)], title="Quadratic errors Warfa")

# err_rwm[,2:8] <- sqrt(err_rwm[,2:8])
# err_mix[,2:8] <- sqrt(err_mix[,2:8])

# graphConvMC_diff2(a[,c(1,3,6,9)],b[,c(1,3,6,9)], title="Quadratic errors Warfa")
# graphConvMC_diff2(err_rwm[-1,c(1,3,6,9)],err_mix[-1,c(1,3,6,9)], title="Quadratic errors Warfa")

# graphConvMC_diff2(a[,c(1,3,6,9)],b[,c(1,3,6,9)])

c <- graphConvMC_diff3(a[,c(1,2,9)],b[,c(1,2,9)])
d <- graphConvMC_diff3(a[,c(1,5,9)],b[,c(1,5,9)])

grid.arrange(c,d, ncol=2)
