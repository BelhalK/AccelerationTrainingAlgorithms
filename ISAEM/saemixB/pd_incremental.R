
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


PD1.saemix<-read.table( "data/PD1.saemix.tab",header=T,na=".")
PD1.saemix <- subset(PD1.saemix, dose!="90")


PD2.saemix<-read.table( "data/PD2.saemix.tab",header=T,na=".")
saemix.data1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))
saemix.data2<-saemixData(name.data=PD2.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
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


saemix.model<-saemixModel(model=modelemax,description="Emax model",type="structural",
psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")



K1 = 500
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2
batchsize<-50

#Weibull
options_pd<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
pd<-data.frame(saemix(saemix.model,saemix.data2,options_pd))
pd<-cbind(iterations,pd)

options_pdincr<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize)
pdincr<-data.frame(saemix_incremental(saemix.model,saemix.data2,options_pdincr))
pdincr<-cbind(iterations,pdincr)

graphConvMC2_saem(pd,pdincr, title="new kernel")


pd$algo <- 'rwm'
pdincr$algo <- 'ISAEM'

pd_scaled <- pd[rep(seq_len(nrow(pd)), each=100/batchsize),]
pd_scaled$iterations = 1:(2*(K1+K2+1))


comparison <- 0
# comparison <- rbind(theo_ref,theo_incremental)
comparison <- rbind(pd_scaled[iterations,],pdincr)

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)

options_pdnew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0,map.range=c(1:5))
pdnew<-data.frame(saemix(saemix.model,saemix.data2,options_pdnew))



final_rwm <- 0
final_ref <- 0
var_rwm <- 0
error_rwm <- 0
final_mix <- 0
final_mix25 <- 0
var_mix <- 0
var_mix25 <- 0
error_mix <- 0
error_mix25 <- 0


lambda_true <- 2
o_lambda_true <- 0.3
beta_true <- 2
o_beta_true <- 0.3
final_mix <- 0
true_param <- c(lambda_true,beta_true,o_lambda_true,o_beta_true)


seed0 = 39546
replicate = 20


for (j in 1:replicate){

     
model2 <- inlineModel("
                      [LONGITUDINAL]
                      input = {e0, emax, e50}

                      
                      EQUATION:
                      Cc = e0+emax*t/(e50+t)
                      
                      DEFINITION:
                      y1 ={distribution=normal, prediction=Cc, sd=a}
                      
                      [INDIVIDUAL]
                      input={e0_pop,o_e0,emax_pop,o_emax,e50_pop,o_e50}
                      
                      DEFINITION:
                      e0  ={distribution=lognormal, prediction=e0_pop,  sd=o_e0}
                      emax   ={distribution=lognormal, prediction=emax_pop,   sd=o_emax}
                      e50  ={distribution=lognormal, prediction=e50_pop,  sd=o_e50}
                      ")

adm  <- list(amount=1, time=seq(0,50,by=50))

p <- c(e0_pop=1, o_e0=0.5,
       emax_pop=10, o_emax=0.3, 
       e50_pop=0.2, o_e50=0.3,  
       a=0.1)
y1 <- list(name='y1', time=seq(1,to=50,by=2))


res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=100, level="individual"),
                 output = y1)



writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_pd.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_pd.csv", header=T, sep=","))
  pd.data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_pd.csv", header=T, sep=",")
  pd.data <- pd.data[pd.data$ytype==2,]

  saemix.data<-saemixData(name.data=pd.data,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))


saemix.model<-saemixModel(model=timetoevent.model,description="time model",type="likelihood",   
  psi0=matrix(c(2,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))
  print(j)

  options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=100)
  theo_ref<-data.frame(saemix_incremental(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  var_rwm <- var_rwm + (theo_ref[,2:5]-true_param)^2
  ML <- theo_ref[,2:5]
  ML[1:(end+1),]<- theo_ref[end+1,2:5]
  error_rwm <- error_rwm + (theo_ref[,2:5]-ML)^2
  theo_ref['individual'] <- j
  final_ref <- rbind(final_ref,theo_ref)


  options.incremental<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize50)
  theo_mix<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental))
  theo_mix <- cbind(iterations, theo_mix)

  var_mix <- var_mix + (theo_mix[,2:5]-true_param)^2
  ML <- theo_mix[,2:5]
  ML[1:(end+1),]<- theo_mix[end+1,2:5]
  error_mix <- error_mix + (theo_mix[,2:5]-ML)^2
  theo_mix['individual'] <- j
  final_mix <- rbind(final_mix,theo_mix)
  
  options.incremental25<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), nb.replacement=batchsize25)
  theo_mix25<-data.frame(saemix_incremental(saemix.model,saemix.data,options.incremental25))
  theo_mix25 <- cbind(iterations, theo_mix25)
  
  var_mix25 <- var_mix25 + (theo_mix25[,2:5]-true_param)^2
  ML <- theo_mix25[,2:5]
  ML[1:(end+1),]<- theo_mix25[end+1,2:5]
  error_mix25 <- error_mix25 + (theo_mix25[,2:5]-ML)^2
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


err_rwm[,2:5] <- error_rwm[,2:5]
err_mix[,2:5] <- error_mix[,2:5]
err_mix25[,2:5] <- error_mix25[,2:5]

err_mix[2,] = err_rwm[2,]=err_mix25[2,]

err_rwm_scaled <- err_rwm[rep(seq_len(nrow(err_rwm)), each=4),]
err_mix_scaled <- err_mix[rep(seq_len(nrow(err_mix)), each=2),]
err_rwm_scaled$iterations = 1:(4*(K1+K2+1))
err_mix_scaled$iterations = 1:(2*(K1+K2+1))


# c <- graphConvMC_se2(err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)],err_rwm_scaled[,c(1,2,8)])
c <- graphConvMC_sec(err_rwm_scaled[2:end,c(1,2,6)],err_mix_scaled[2:end,c(1,2,6)],err_mix25[2:end,c(1,2,6)])
d <- graphConvMC_sed(err_rwm_scaled[2:end,c(1,4,6)],err_mix_scaled[2:end,c(1,4,6)],err_mix25[2:end,c(1,4,6)])

grid.arrange(c,d, ncol=2)


