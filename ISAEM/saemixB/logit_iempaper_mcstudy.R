setwd("/Users/karimimohammedbelhal/Desktop/CSDA_code/Dir")
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
  source('main_time.R')
  source('main_estep_time.R')
  source('main_mstep_time.R') 
  source('func_aux_time.R') 
  source('SaemixObject_time.R') 
  source('main_initialiseMainAlgo_time.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat_incremental")
source('main_cat2.R')
source('main_estep_cat2.R')
source('main_cat_incremental.R')
source('estep_cat_incremental.R')
source('main_mstep_cat.R') 
source('func_aux_cat.R') 
source('SaemixObject_cat.R') 
source('main_initialiseMainAlgo_cat.R') 
source("mixtureFunctions.R")

source("/Users/karimimohammedbelhal/Desktop/papers/iem_code/imcem_saemix/plots_se.R")


library("rJava")
library("rCMA")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
iter_mcmc = 200


# cat_data.saemix<-read.table("data/categorical1_data.txt",header=T,na=".")
# cat_data.saemix <- cat_data.saemix[1:692,]
# saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("ID"),name.response=c("Y"),name.predictors=c("Y"), name.X=c("TIME"))



# cat_data.saemix<-read.table("data/categorical1_data.txt",header=T,na=".")
# cat_data.saemix<-read.table("data/categorical1_data_less2.txt",header=T,na=".")
# cat_data.saemix<-read.table("data/categorical1_data_less2.txt",header=T,na=".")
# cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat1.csv", header=T, sep=",")
cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat2.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("y"), name.X=c("time"))
# saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("ID"),name.response=c("Y"),name.predictors=c("Y"), name.X=c("TIME"))



cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]

th1 <- psi[id,1]
th2 <- psi[id,2]
th3 <- psi[id,3]

P0 <- 1/(1+exp(-th1))
Pcum1 <- 1/(1+exp(-th1-th2))
Pcum2 <- 1/(1+exp(-th1-th2-th3))

P1 <- Pcum1 - P0
P2 <- Pcum2 - Pcum1
P3 <- 1 - Pcum2

P.obs = (level==0)*P0+(level==1)*P1+(level==2)*P2+(level==3)*P3

return(P.obs)
}


saemix.model<-saemixModel(model=cat_data.model,description="cat model",   
  psi0=matrix(c(2,1,2),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))), 
  transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),omega.init=matrix(c(2,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")


saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0))
saemix.foce<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc))


# post_rwm<-saemix_post_cat(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
# post_foce<-saemix_post_cat(saemix.model,saemix.data,saemix.foce)$post_newkernel


K1 = 500
K2 = 100

iteration = 1:(K1+K2+1)
gd_step = 0.01
end = K1+K2
seed0 = 444

replicate = 3

final_rwm <- 0
final_incremental <- 0
final_incremental25 <- 0
for (m in 1:replicate){
  print(m)
  print(m)
  # l = list(c(1,5,1,0,0,0),c(0.8,12,0.8,0,0,0),c(1.2,3,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  l = list(c(2,1,2),c(2,1,5),c(2,1,3))
  saemix.model<-saemixModel(model=cat_data.model,description="cat model",   
  psi0=matrix(l[[m]],ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))), 
  transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),omega.init=matrix(c(3/m,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")


  options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=100)
  theo_ref<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options))
  theo_ref <- cbind(iteration, theo_ref)
  theo_ref['individual'] <- m
  theo_ref_scaled <- theo_ref[rep(seq_len(nrow(theo_ref)), each=4),]
  theo_ref_scaled$iteration = 1:(4*(K1+K2+1))
  final_rwm <- rbind(final_rwm,theo_ref_scaled[iteration,])

  options.incremental<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=50)
  theo_mix<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental))
  theo_mix <- cbind(iteration, theo_mix)
  theo_mix['individual'] <- m
  theo_mix_scaled <- theo_mix[rep(seq_len(nrow(theo_mix)), each=2),]
  theo_mix_scaled$iteration = 1:(2*(K1+K2+1))
  final_incremental <- rbind(final_incremental,theo_mix_scaled[iteration,])

   options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=25)
  theo_mix25<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental25))
  theo_mix25 <- cbind(iteration, theo_mix25)
  theo_mix25['individual'] <- m
  final_incremental25 <- rbind(final_incremental25,theo_mix25[,])
}

graphConvMC_diff3(final_incremental,final_incremental, final_incremental,title="Diff intial param Warfa")

graphConvMC_diff2(final_rwm[,c(1,3,9)],final_incremental[,c(1,3,9)], title="Diff intial param Warfa")

graphConvMC_diff2(final_rwm[,c(1,3,6,9)],final_incremental[,c(1,3,6,9)], title="Diff intial param Warfa")

graphConvMC_diff(final_rwm[,c(1,2,9)],final_incremental[,c(1,2,9)], title="Diff intial param Warfa")

graphConvMC_diff2(final_rwm[,c(1,3,6,9)],final_incremental[,c(1,3,6,9)])




# a <- graphConvMC_diff4(final_rwm[,c(1,4,8)],final_incremental[,c(1,4,8)])
# b <- graphConvMC_diff3(final_rwm[,c(1,7,8)],final_incremental[,c(1,7,8)])

# grid.arrange(a,b, ncol=2)
colnames(final_rwm)[4] <- "z3"
colnames(final_incremental)[4] <- "z3"
colnames(final_incremental25)[4] <- "z3"
a <- graphConvMC_diff4(final_rwm[,c(1,4,8)],final_incremental[,c(1,4,8)],final_incremental25[,c(1,4,8)])
b <- graphConvMC_diff3(final_rwm[,c(1,5,8)],final_incremental[,c(1,5,8)],final_incremental25[,c(1,5,8)])

grid.arrange(a,b, ncol=2)




#MC STUDY


final_rwm <- 0
final_ref <- 0
var_rwm <- 0
error_rwm <- 0

th1 <- 1
th2 <- 3
th3 <- 1
o_th1 <- 1
o_th2 <- 0.5
o_th3 <- 0.1
final_mix <- 0
final_mix25 <- 0
true_param <- c(th1,th2,th3,o_th1,o_th2,o_th3)
var_mix <- 0
var_mix25 <- 0
error_mix <- 0
error_mix25 <- 0


seed0 = 39546
replicate = 30

for (j in 1:replicate){

     model2 <- inlineModel("

                  [LONGITUDINAL]
                  input = {th1, th2, th3}

                  EQUATION:
                  lgp0 = th1
                  lgp1 = lgp0 + th2
                  lgp2 = lgp1 + th3

                  DEFINITION:
                  level = { type = categorical,  categories = {0, 1, 2, 3},
                  logit(P(level<=0)) = th1
                  logit(P(level<=1)) = th1 + th2
                  logit(P(level<=2)) = th1 + th2 + th3
                  }

                  [INDIVIDUAL]
                  input={th1_pop, o_th1,th2_pop, o_th2,th3_pop, o_th3}
                          

                  DEFINITION:
                  th1  ={distribution=normal, prediction=th1_pop,  sd=o_th1}
                  th2  ={distribution=lognormal, prediction=th2_pop,  sd=o_th2}
                  th3  ={distribution=lognormal, prediction=th3_pop,  sd=o_th3}
                          
                          ")

    p <- c(th1_pop=2, o_th1=0.8,
           th2_pop=3, o_th2=0.5, 
           th3_pop=1, o_th3=0.7)



    y1 <- list(name='level', time=seq(1,to=100,by=10))



    res2a2 <- simulx(model = model2,
                     parameter = p,
                     group = list(size=100, level="individual"),
                     output = y1)


    writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat_mcstudy.csv")


  cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat_mcstudy.csv", header=T, sep=",")
  saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("y"), name.X=c("time"))
  saemix.model<-saemixModel(model=cat_data.model,description="cat model",   
  psi0=matrix(c(1,1,1),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))), 
  transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3, 
  byrow=TRUE),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")
  print(j)
  options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=100)
  theo_ref<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options))
  theo_ref <- cbind(iteration, theo_ref)
  var_rwm <- var_rwm + (theo_ref[,2:6]-true_param)^2
  ML <- theo_ref[,2:6]
  ML[1:(end+1),]<- theo_ref[end+1,2:6]
  error_rwm <- error_rwm + (theo_ref[,2:6]-ML)^2
  theo_ref['individual'] <- j
  final_ref <- rbind(final_ref,theo_ref)


  options.incremental<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=50)
  theo_mix<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental))
  theo_mix <- cbind(iteration, theo_mix)

  var_mix <- var_mix + (theo_mix[,2:6]-true_param)^2
  ML <- theo_mix[,2:6]
  ML[1:(end+1),]<- theo_mix[end+1,2:6]
  error_mix <- error_mix + (theo_mix[,2:6]-ML)^2
  theo_mix['individual'] <- j
  final_mix <- rbind(final_mix,theo_mix)
  
  options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=25)
  theo_mix25<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental25))
  theo_mix25 <- cbind(iteration, theo_mix25)
  
  var_mix25 <- var_mix25 + (theo_mix25[,2:6]-true_param)^2
  ML <- theo_mix25[,2:6]
  ML[1:(end+1),]<- theo_mix25[end+1,2:6]
  error_mix25 <- error_mix25 + (theo_mix25[,2:6]-ML)^2
  theo_mix25['individual'] <- j
  final_mix25 <- rbind(final_mix25,theo_mix25)
}

graphConvMC_diff(final_ref,final_ref,final_ref)
graphConvMC_diff(final_ref,final_mix,final_mix25)

error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix
error_mix25 <- 1/replicate*error_mix25

error_rwm <- cbind(iterations, error_rwm)
error_mix <- cbind(iterations, error_mix)
error_mix25 <- cbind(iterations, error_mix25)

err_mix<- theo_ref
err_rwm<- theo_ref
err_mix25<- theo_ref


err_rwm[,2:6] <- error_rwm[,1:5]
err_mix[,2:6] <- error_mix[,1:5]
err_mix25[,2:6] <- error_mix25[,1:5]

err_mix[2,] = err_rwm[2,]=err_mix25[2,]

err_rwm_scaled <- err_rwm[rep(seq_len(nrow(err_rwm)), each=4),]
err_mix_scaled <- err_mix[rep(seq_len(nrow(err_mix)), each=2),]
err_rwm_scaled$iteration = 1:(4*(K1+K2+1))
err_mix_scaled$iteration = 1:(2*(K1+K2+1))

c <- graphConvMC_se2(err_rwm_scaled[-1,c(1,2,8)],err_mix_scaled[-1,c(1,2,8)],err_mix25[-1,c(1,2,8)])
d <- graphConvMC_se2(err_rwm_scaled[-1,c(1,3,8)],err_mix_scaled[-1,c(1,3,8)],err_mix25[-1,c(1,3,8)])

grid.arrange(c,d, ncol=2)

# c <- graphConvMC_se2(err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)])
c <- graphConvMC_sec(err_rwm_scaled[2:600,c(1,3,7)],err_mix_scaled[2:600,c(1,3,7)],err_mix25[2:600,c(1,3,7)])
d <- graphConvMC_sed(err_rwm_scaled[2:600,c(1,6,7)],err_mix_scaled[2:600,c(1,6,7)],err_mix25[2:600,c(1,6,7)])

grid.arrange(c,d, ncol=2)

