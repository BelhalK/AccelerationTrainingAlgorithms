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
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/")
source('plots.R') 

library('rCMA')
###WARFA
warfa_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/data/warfarin_data.txt", header=T)
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
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

saemix.model_warfanovar<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE))


K1 = 200
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2

#With var 
options_warfa_without<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_without<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa_without))
warfa_without <- cbind(iterations, warfa_without)

options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=0)
warfa_newkernel<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_newkernel))
warfa_newkernel <- cbind(iterations, warfa_newkernel)

#No var
##### Optim (fmin search)
options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
warfa_optim<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
warfa_optim <- cbind(iterations, warfa_optim)

##### AV
options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=1)
warfa_withav<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
warfa_withav <- cbind(iterations, warfa_withav)


##### AV and new kernel
options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5), av=1)
warfa_newkernelav<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_newkernel))
warfa_newkernelav <- cbind(iterations, warfa_newkernelav)

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


options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0,map.range=c(0))
warfa_bayes<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
warfa_bayes <- cbind(iterations, warfa_bayes)



replicate = 3
seed0 = 395246

#RWM
final_optim <- 0
final_av <- 0
final_avnew <- 0
final_bayes <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(1,7,1,0,0,0),c(0.8,7.2,0.8,0,0,0),c(1.2,6.8,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  
  saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

  saemix.model_warfanovar<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
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
  options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0)
  warfa_optim<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
  warfa_optim <- cbind(iterations, warfa_optim)
  warfa_optim['individual'] <- m
  final_optim <- rbind(final_optim,warfa_optim)
  ##### AV
  options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, map.range=c(0),av=1)
  warfa_withav<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
  warfa_withav <- cbind(iterations, warfa_withav)
  warfa_withav['individual'] <- m
  final_av <- rbind(final_av,warfa_withav)

  ##### AV and new kernel
  options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0,map.range=c(1:3), av=1)
  warfa_newkernelav<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_newkernel))
  warfa_newkernelav <- cbind(iterations, warfa_newkernelav)
  warfa_newkernelav['individual'] <- m
  final_avnew <- rbind(final_avnew,warfa_newkernelav)

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


  options_warfa_with<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,2),nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, av=0,map.range=c(0))
  warfa_bayes<-data.frame(saemix(saemix.model_warfanovar,saemix.data_warfa,options_warfa_with))
  warfa_bayes <- cbind(iterations, warfa_bayes)
  warfa_bayes['individual'] <- m
  final_bayes <- rbind(final_bayes,warfa_bayes)
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

graphConvMC_diff(final_av,final_avnew, title="")

#black: optim
#blue: av
#red: av newkernel
#green: bayes