setwd("/Users/karimimohammedbelhal/Desktop/monolix/SEAMIX ODE/Dir2")
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_cov.R') 
  source('func_distcond.R') 
  source('func_FIM.R') 
  source('func_ggplot2.R') 
  source('func_plots.R') 
  source('func_simulations.R') 
  source('ggplot2_global.R') 
  source('global.R') 
  source('main.R') 
  source('main_estep.R') 
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('plots_ggplot2.R') 
  source('saemix-package.R') 
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
source('MlxCore.R') 
source('APIManager.R') 
source('APITools.R') 
source('CovariateModels.R') 
source('displayTools.R') 
source('IndividualModel.R') 
source('MlxEnvironment.R') 
source('ObservationModels.R') 
source('PopulationParameters.R') 
source('ProjectManagement.R') 
source('Results.R') 
source('Scenario.R')   
source('Settings.R')   
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/cmaes/")
source("mixtureFunctions.R")


library("rJava")
library("rCMA")
library("mlxR")
library(sgd)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)


#####################################################################################
# Theophylline

# Data - changing gender to M/F
theo.saemix<-read.table("theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")



# theo.saemix_less <- theo.saemix[1:120,]
theo.saemix_less <- theo.saemix
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
saemix.data<-saemixData(name.data=theo.saemix_less,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")


modelMM<- function(psi, id, xidep){
  Model <- inlineModel("
                       [LONGITUDINAL]
                       input = {kf,kr,kcat}
                       EQUATION:
                       E_0=200
                       S_0=500
                       ES_0=0
                       f_0=0
                       t_0=0
                       ddt_E=-kf*E*S+kr*ES+kcat*ES
                       ddt_S=-kf*E*S+kr*ES
                       ddt_ES=kf*E*S-kr*ES-kcat*ES
                       ddt_f=kcat*ES
                             ")
  time <- rep(0,length(id))
  p <- data.frame(id=1:id[length(id)],kf=psi[,1],kr=psi[,2],kcat=psi[,3])
  time <- xidep[,1]
  design.f <- as.data.frame(cbind(id,time))
  f <- list(name="f", time=design.f)
  resList <- list("Model"=Model,"p"=p,"f"=f) 
  return(resList)
  
}


saemix.model<-saemixModel(model=modelMM,
                          description="Michaelis Menten", 
                          psi0=matrix(c(0.001,0.0001,0.1),ncol=3, byrow=TRUE,
                          dimnames=list(NULL, c("kf","kr","kcat"))),transform.par=c(1,1,1),
                          fixed.estim=c(1,1,1),
                          covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
                          omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")

saemix.fit<-saemix(saemix.model,saemix.data,list(seed=632545,directory="MichaelisMenten",
                                                 save=FALSE,save.graphs=FALSE,map=FALSE,fim=FALSE,ll.is=FALSE))
