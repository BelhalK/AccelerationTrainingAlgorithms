setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/saemixrandomvariable")
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
  source('SaemixObject.R') 
  source('zzz.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/")
source('plots.R') 

library('rCMA')
###logit


logit_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/data/Data_logistic_1D_simule.txt", header=T)
saemix.data_logit<-saemixData(name.data=logit_data,header=TRUE,sep=" ",na=NA, name.group=c("group"),
  name.predictors=c("X"),name.response=c("Y"), name.X="X")

model1cpt<-function(psi,id,xidep) { 
  tim<-xidep[,1]  
  p0<-psi[id,1]
  alpha<-psi[id,2]
  tau<-psi[id,3]
  ypred<-1/(1+((1/p0)-1)*exp(-alpha*(tim-tau)/(p0*(1-p0))))
  return(ypred)
}


saemix.model_logit<-saemixModel(model=model1cpt,description="logitrin",type="structural"
  ,psi0=matrix(c(0.6,0.05,60),ncol=3,byrow=TRUE, dimnames=list(NULL, c("p0","alpha","tau"))),
  transform.par=c(3,1,0),omega.init=matrix(c(1,0,0,0,1,0,0,0,100),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

saemix.model_logitnovar<-saemixModel(model=model1cpt,description="logitrin",type="structural"
  ,psi0=matrix(c(0.6,0.05,60),ncol=3,byrow=TRUE, dimnames=list(NULL, c("p0","alpha","tau"))),
  transform.par=c(3,1,0),omega.init=matrix(c(1,0,0,0,1,0,0,0,100),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(0,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


K1 = 200
K2 = 10
iterations = 1:(K1+K2+1)
end = K1+K2


#With var no sa
options_logit_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
logit_without<-data.frame(saemix(saemix.model_logit,saemix.data_logit,options_logit_without))
logit_without <- cbind(iterations, logit_without)

graphConvMC_twokernels(logit_without,logit_without)


options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
logit_newkernel<-data.frame(saemix(saemix.model_logit,saemix.data_logit,options_newkernel))
logit_newkernel <- cbind(iterations, logit_newkernel)


#no var no sa
options_logit_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
logit_without<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_logit_without))
logit_without <- cbind(iterations, logit_without)

graphConvMC_twokernels(logit_without,logit_without)


options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
logit_newkernel<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_newkernel))
logit_newkernel <- cbind(iterations, logit_newkernel)


#No var no sa but randomvariable
options_logit_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0,map.range=c(0))
logit_withnosa<-data.frame(saemix(saemix.model_logitnovar,saemix.data_logit,options_logit_with))
logit_withnosa <- cbind(iterations, logit_withnosa)
graphConvMC_twokernels(logit_without[,1:6],logit_withnosa)

options_logit_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0,map.range=c(0))
logit_withnosa<-data.frame(saemix(saemix.model_logitnovar_init,saemix.data_logit,options_logit_with))
logit_withnosa <- cbind(iterations, logit_withnosa)



replicate = 6
seed0 = 395246

#RWM
final_optim <- 0
final_av <- 0
final_avnew <- 0
final_bayes <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(0.5,0.05,60),c(0.55,0.07,50),c(0.6,0.1,40),c(0.15,0.05,60),c(0.20,0.05,60),c(0.25,0.05,60))
  
  saemix.model_logit<-saemixModel(model=model1cpt,description="logitrin",type="structural"
  ,psi0=matrix(l[[m]],,ncol=3,byrow=TRUE, dimnames=list(NULL, c("p0","alpha","tau"))),
  transform.par=c(3,1,0),omega.init=matrix(c(1,0,0,0,1,0,0,0,100),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

saemix.model_logisticnovar<-saemixModel(model=model1cpt,description="logitrin",type="structural"
  ,psi0=matrix(l[[m]],,ncol=3,byrow=TRUE, dimnames=list(NULL, c("p0","alpha","tau"))),
  transform.par=c(3,1,0),omega.init=matrix(c(1,0,0,0,1,0,0,0,100),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(0,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


  setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/saemixnonvariability")
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
  source('SaemixObject.R') 
  source('zzz.R') 

  #No var
  ##### Optim (fmin search)
  options_logistic_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
  logistic_optim<-data.frame(saemix(saemix.model_logisticnovar,saemix.data_logit,options_logistic_with))
  logistic_optim <- cbind(iterations, logistic_optim)
  logistic_optim['individual'] <- m
  final_optim <- rbind(final_optim,logistic_optim)
  ##### AV
  options_logistic_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, map.range=c(0),av=1)
  logistic_withav<-data.frame(saemix(saemix.model_logisticnovar,saemix.data_logit,options_logistic_with))
  logistic_withav <- cbind(iterations, logistic_withav)
  logistic_withav['individual'] <- m
  final_av <- rbind(final_av,logistic_withav)

  ##### AV and new kernel
  options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6),nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1), av=1)
  logistic_newkernelav<-data.frame(saemix(saemix.model_logisticnovar,saemix.data_logit,options_newkernel))
  logistic_newkernelav <- cbind(iterations, logistic_newkernelav)
  logistic_newkernelav['individual'] <- m
  final_avnew <- rbind(final_avnew,logistic_newkernelav)

  ##### pseudo bayesian
  setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/saemixrandomvariable")
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
    source('SaemixObject.R') 
    source('zzz.R') 


  options_logistic_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0,map.range=c(0))
  logistic_bayes<-data.frame(saemix(saemix.model_logisticnovar,saemix.data_logit,options_logistic_with))
  logistic_bayes <- cbind(iterations, logistic_bayes)
  logistic_bayes['individual'] <- m
  final_bayes <- rbind(final_bayes,logistic_bayes)
}


graphConvMC_diff <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+
      xlab("iteration")+scale_x_log10()+ ylab(names(df[j]))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}


graphConvMC_diff4 <- function(df,df2,df3,df4, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  df3$individual <- as.factor(df3$individual)
  df4$individual <- as.factor(df4$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red")+geom_line(aes_string(df4[,1],df4[,j],by=df4[,ncol(df4)]),colour="green")+
      xlab("iteration") +scale_x_log10()+ ylab(names(df[j]))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}

graphConvMC_diff4(final_optim,final_av,final_avnew,final_bayes, title="")

#black: optim
#blue: av
#red: av newkernel
#green: bayes

graphConvMC_diff(final_optim,final_bayes, title="")

graphConvMC_diff(subset(final_optim, individual==c(4,5,6)),subset(final_bayes, individual==c(4,5,6)), title="")
graphConvMC_diff4(subset(final_optim, individual==c(1,2,3)),subset(final_bayes, individual==c(1,2,3)),subset(final_optim, individual==c(4,5,6)),subset(final_bayes, individual==c(4,5,6)), title="")
