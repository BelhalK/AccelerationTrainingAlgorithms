setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Annealing - Main/saemixConnectorsMonolix/R")
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
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/Annealing - Main/saemixConnectorsMonolix")

library("mlxR")
library(MlxConnectors)
initializeMlxConnectors(software = "monolix")

################################################################ MonolixProject ####################################################################################################################################

project.file <- "mlxProjects/hcv/hcv_project.mlxtran"
loadProject(project.file)

getEstimatedPopulationParameters()
computePredictions(getEstimatedIndividualParameters()$saem)
computePredictions(getEstimatedIndividualParameters()$saem, individualIds = c(10,20))

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  s<-psi[id,1]
  d<-psi[id,2]
  beta<-psi[id,3]
  delta<-psi[id,4]
  p<-psi[id,5]
  c<-psi[id,6]
  eta<-psi[id,7]
  epsilon<-psi[id,8]
  ypred<-1
  return(ypred)
}

hcv_data <- readDatamlx(project = project.file)
treat <- hcv_data$treatment[,c(3)]
# hcv.saemix <- cbind(hcv_data$Y,treat)
hcv.saemix <- hcv_data$Y

saemix.data<-saemixData(name.data=hcv.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("time"),name.response=c("Y"), name.X="time")

# saemix.model<-saemixModel(model=model1cpt,description="hcv",type="structural"
#   ,psi0=matrix(c(1000,1,0.00005,0.05,20,5,0.9,0.7),ncol=8,byrow=TRUE, dimnames=list(NULL, c("s","d","beta","delta","p","c","eta","epsilon"))),fixed.estim=c(1,1,1,1,1,1,1,1),
#   transform.par=c(1,1,1,1,1,1,1,1),omega.init=matrix(diag(8),ncol=8,byrow=TRUE),covariance.model=matrix(diag(8),ncol=8,byrow=TRUE))

# K1 = 200
# K2 = 100
# iterations = 1:(K1+K2+1)
# end = K1+K2

# options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,nb.chains=1,monolix=TRUE)
# hcv<-saemix(saemix.model,saemix.data,options)



K1 = 100
K2 = 50
iterations = 1:(K1+K2+1)
end = K1+K2
nchains = 1
runtime = 20
replicate = 3
seed0 = 395246

#RWM
final_optim <- 0
final_av <- 0
final_avnew <- 0
final_bayes <- 0

final_annealing <- 0

cov.model = matrix(0,nrow=8,ncol=8,byrow=TRUE)
cov.model[1,1] <- 1
cov.model[2,2] <- 1

for (m in 1:replicate){
  print(m)
  l = list(c(1000,1,0.00005,0.05,20,5,0.9,0.7),
          c(1000,1,0.00005,0.05,20,5,0.9,0.7),
          c(1000,1,0.00005,0.05,20,5,0.9,0.7))
  
  saemix.model_hcvnovar<-saemixModel(model=model1cpt,description="hcv",type="structural"
  ,psi0=matrix(l[[m]],ncol=8,byrow=TRUE, 
    dimnames=list(NULL, c("s","d","beta","delta","p","c","eta","epsilon"))),
  fixed.estim=c(1,1,1,1,1,1,1,1),
  transform.par=c(1,1,1,1,1,1,1,1),omega.init=matrix(diag(8),ncol=8,byrow=TRUE),
  covariance.model=cov.model)


  #No var
  ##### Optim (fmin search)
  options_hcv_with<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
   duration = runtime,an=FALSE,coeff=1,nbiter.burn =0, av=0)
  hcv_optim<-data.frame(saemix(saemix.model_hcvnovar,saemix.data,options_hcv_with))
  hcv_optim <- cbind(iterations, hcv_optim)
  hcv_optim['individual'] <- m
  final_optim <- rbind(final_optim,hcv_optim)
  
  #### AV
  options_hcv_with<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=K1,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
   duration = runtime,an=FALSE,coeff=1,nbiter.burn =0, av=0)
  hcv_withav<-data.frame(saemix(saemix.model_hcvnovar,saemix.data,options_hcv_with))
  hcv_withav <- cbind(iterations, hcv_withav)
  hcv_withav['individual'] <- m
  final_av <- rbind(final_av,hcv_withav)

  # ##### AV and new kernel
  # options_newkernel<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6,0),
  # nbiter.sa=K1/2,nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0,map.range=c(1), av=1)
  # hcv_newkernelav<-data.frame(saemix(saemix.model_hcvnovar,saemix.data,options_newkernel))
  # hcv_newkernelav <- cbind(iterations, hcv_newkernelav)
  # hcv_newkernelav['individual'] <- m
  # final_avnew <- rbind(final_avnew,hcv_newkernelav)


  # ##### Annealing MCMC
  options_an<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
   duration = runtime, an=TRUE,coeff=0.005,
  map.range=c(0), av=1)
  hcv_an<-data.frame(saemix(saemix.model_hcvnovar,saemix.data,options_an))
  hcv_an <- cbind(iterations, hcv_an)
  hcv_an['individual'] <- m
  final_annealing<- rbind(final_annealing,hcv_an)

  ##### pseudo bayesian

  options_hcv_with<-list(seed=39546,map=F,fim=F,ll.is=F,
  nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=FALSE,nbiter.burn =0,nb.chains=nchains,monolix=TRUE,
   duration = runtime,an=FALSE,coeff=1, 
    nbiter.burn =0, av=0,map.range=c(0))
  hcv_bayes<-data.frame(saemix(saemix.model_hcvnovar,saemix.data,options_hcv_with))
  hcv_bayes <- cbind(iterations, hcv_bayes)
  hcv_bayes['individual'] <- m
  final_bayes <- rbind(final_bayes,hcv_bayes)

}


graphConvMC_diff4(final_optim,final_av,final_annealing,final_bayes, title="")
#black: optim
#blue: av
#red: annealing
#green: bayes



