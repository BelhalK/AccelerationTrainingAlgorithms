library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
library(data.table)
library(rstan)

# save.image("RData/mcmc_tobramycin.RData")


source('R/aaa_generics.R') 
source('R/compute_LL.R') 
source('R/func_aux.R') 
source('R/func_distcond.R') 
source('R/func_FIM.R')
source('R/func_plots.R') 
source('R/func_simulations.R') 
source('R/estep_mcmc.R')
source('R/indiv_VI.R')
source('R/variationalinferencelinear.R')
source('R/main.R')
source('R/main_estep.R')
source('R/mcmc_final.R')
source('R/mcmc_indiv.R')
source('R/main_initialiseMainAlgo.R') 
source('R/main_mstep.R') 
source('R/check_linearvslaplace.R')
source('R/SaemixData.R')
source('R/SaemixModel.R') 
source('R/SaemixRes.R') 
# source('R/SaemixRes_c.R') 
source('R/SaemixObject.R') 
source('R/zzz.R')
source('R/graphplot.R')




# bolus_data <- read.table("data/bolus_data.txt", header=T)
# saemix.data<-saemixData(name.data=bolus_data,header=TRUE,sep=" ",na=NA, name.group=c("ID"),
#   name.predictors=c("amt","time"),name.response=c("y"), name.X="time")

bolus_data <- read.table("data/tobramycin.txt", header=T)
bolus_data <- bolus_data[bolus_data$EVID==0,]
saemix.data<-saemixData(name.data=bolus_data,header=TRUE,sep=" ",na=NA, name.group=c("ID"),
  name.predictors=c("DOSE","TIME"),name.response=c("CP"), name.X="TIME")



model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  V<-psi[id,1]
  k<-psi[id,2]
  CL<-k*V
  # ypred<-dose/(V*k)*exp(-k*tim)
  ypred<-1/(V*k)*exp(-k*tim)
  return(ypred)
}


##RUNS
K1 = 400
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

# Default model, no covariate
saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural",
  ,psi0=matrix(c(70,1),ncol=2,byrow=TRUE, dimnames=list(NULL, c("V","k"))),
  transform.par=c(1,1),omega.init=matrix(c(1,0,0,1),ncol=2,byrow=TRUE),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))


i = 2

L_mcmc=20000
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=i)
ref<-mcmc.indiv(saemix.model_warfa,saemix.data,options_warfa)$eta


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=i)
new<-mcmc.indiv(saemix.model_warfa,saemix.data,options_warfanew)$eta

new.student<-mcmc.indiv(saemix.model_warfa,saemix.data,options_warfanew)$eta

#MALA

options.mala<-list(seed=39546,map=F,fim=F,ll.is=F, av=0, sigma.val=0.002
  ,gamma.val=0.1,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = i)
mala<-mcmc.indiv(saemix.model_warfa,saemix.data,options.mala)$eta


#RSTAN VB

model <- 'data {
          int<lower=0> N;// Number of observations
          vector[N] time; //predictor
          real dose; //predictor
          vector[N] concentration;  //response
          
          real beta1_pop;
          real beta2_pop;
          real beta3_pop;
          real<lower=0> omega_beta1;
          real<lower=0> omega_beta2;
          real<lower=0> omega_beta3;
          real<lower=0>  pres;
        }
        parameters {
          vector<lower=0>[3] beta;
        }
        model {
          //Priors
          beta[1] ~ lognormal( beta1_pop , omega_beta1);
          beta[2] ~ lognormal( beta2_pop , omega_beta2);
          beta[3] ~ lognormal( beta3_pop , omega_beta3);

          concentration ~ normal(dose*beta[1]/(beta[2]*(beta[1]-beta[3]))*(exp(-beta[3]*time)-exp(-beta[1]*time)), pres);
        }'

modelstan <- stan_model(model_name = "warfarin",model_code = model)

#NUTS using rstan
L_mcmc <- 100000
i <- 10
options.vi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
vi<-mcmc.indiv(saemix.model_warfa,saemix.data,options.vi)$eta


#ADVI for VI post outputs
#Calculate mu and gamma of ELBO optimization
variational.post.options<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nb.chains=1,
 nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),
  modelstan = modelstan, indiv.index = i)

variational.post<-indiv.variational.inference(saemix.model_warfa,saemix.data,variational.post.options)
mu.vi <- variational.post$mu
Gamma.vi <- variational.post$Gamma
etamap <- variational.post$map
Gammamap <- variational.post$Gammamap

# #using the output of ADVI (drawn from candidate KL posterior)
# test <- etamap
# # test[i,] <- etamap[i,] +0.01
# test[i,] <- mu.vi
eta.vi <- etamap
Gammavi <- Gammamap
# eta.vi[i,] <- mu.vi
# Gammavi[[i]] <- Gamma.vi
options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=eta.vi,Gamma = Gammavi,
        nbiter.mcmc = c(0,0,0,0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),
        nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index = i)
advi<-mcmc.indiv(saemix.model_warfa,saemix.data,options_warfavi)$eta

abs(mu.vi - etamap[i,])/abs(etamap[i,])
norm(Gamma.vi - Gammamap[[i]])/norm(Gammamap[[i]])
abs(Gamma.vi - Gammamap[[i]])/abs(Gammamap[[i]])
(Gamma.vi - Gammamap[[i]])/(Gammamap[[i]])


i <- 3
#Autocorrelation
rwm.obj <- as.mcmc(ref)
refac <- autocorr.plot(rwm.obj[,1]) 

new.obj <- as.mcmc(new)
newac <- autocorr.plot(new.obj[,1]) 

newstudent.obj <- as.mcmc(new.student)
newstudentac <- autocorr.plot(newstudent.obj[,1]) 

mala.obj <- as.mcmc(mala)
malaac <- autocorr.plot(mala.obj[,1]) 

vi.obj <- as.mcmc(vi)
viac <- autocorr.plot(vi.obj[,1]) 


advi.obj <- as.mcmc(advi)
adviac <- autocorr.plot(advi.obj[,1]) 

#Autocorrelation
par(mfrow=c(1,5))
acf(ref[,1], main="RWM")
acf(new[,1], main="IMH (Gaussian)")
acf(new.student[,1], main="IMH (Student)")
acf(mala[,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(advi[,1], main="ADVI")

par(mfrow=c(1,5))
acf(ref[,2], main="RWM")
acf(new[,2], main="IMH")
acf(mala[,2], main="MALA")
acf(vi[,2], main="NUTS")
acf(advi[,2], main="ADVI")


par(mfrow=c(1,5))
acf(ref[,3], main="RWM")
acf(new[,3], main="IMH")
acf(mala[,3], main="MALA")
acf(vi[,3], main="NUTS")
acf(advi[,3], main="ADVI")


#MSJD
mssd(ref[,1])
mssd(new[,1])
mssd(new.student[,1])
mssd(mala[,1])
mssd(advi[,1])
mssd(vi[,1])


mssd(refaverage[,1])
mssd(refaverage[,1])
mssd(refaverage[,1])
mssd(refaverage[,1])
mssd(refaverage[,1])


# #Autocorrelation
# rwm.obj <- as.mcmc(states.ref[[10]])
# autocorr.plot(rwm.obj[,1]) + title("RWM Autocorrelation")

# new.obj <- as.mcmc(states.new[[10]])
# autocorr.plot(new.obj[,1]) + title("Laplace Autocorrelation")

# mala.obj <- as.mcmc(states.mala[[10]])
# autocorr.plot(mala.obj[,1]) + title("MALA Autocorrelation")

# nuts.obj <- as.mcmc(states.nuts[[10]])
# autocorr.plot(nuts.obj[,1]) + title("NUTS Autocorrelation")





L_mcmc=20000
i <- 10
nchains <- 50
#REF
listofrefchains <- list(ref,ref)
for (m in 1:nchains){
  options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(2,2,2,0,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=i)
  ref<-mcmc.indiv(saemix.model_warfa,saemix.data,options_warfa)$eta
  listofrefchains[[m]] <- ref
}


#new
listofnewchains <- list(ref,ref)
for (m in 1:nchains){
  print(m)
options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=i)
new<-mcmc.indiv(saemix.model_warfa,saemix.data,options_warfanew)$eta
listofnewchains[[m]] <- new
}


#new
listofnewstudentchains <- list(ref,ref)
for (m in 1:nchains){
  print(m)
options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0,0,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,
  displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index=i)
new.student<-mcmc.indiv(saemix.model_warfa,saemix.data,options_warfanew)$eta
listofnewstudentchains[[m]] <- new.student
}


#MALA
listofmalachains <- list(ref,ref)
for (m in 1:nchains){
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F, av=0, sigma.val=0.002
  ,gamma.val=0.1,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,6,0,0),nb.chains=1
  , nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0
  , map.range=c(0), indiv.index = i)
mala<-mcmc.indiv(saemix.model_warfa,saemix.data,options.mala)$eta
listofmalachains[[m]] <- mala
}


#NUTS
listofnutschains <- list(ref,ref)
# listofnutschains <- 0
for (m in 1:nchains){
options.vi<-list(seed=39546*m,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,
  nbiter.mcmc = c(0,0,0,0,0,1,0),nb.chains=1, nbiter.saemix = c(K1,K2),
  nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), 
  modelstan = modelstan, indiv.index = i)
vi<-mcmc.indiv(saemix.model_warfa,saemix.data,options.vi)$eta
  # listofnutschains <- listofnutschains + nuts
listofnutschains[[m]] <- vi
}



#ADVI
listofadvichains <- list(ref,ref)
for (m in 1:nchains){
options_warfavi<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc, mu=eta.vi,Gamma = Gammavi,
        nbiter.mcmc = c(0,0,0,0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),
        nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0), indiv.index = i)
advi<-mcmc.indiv(saemix.model_warfa,saemix.data,options_warfavi)$eta
listofadvichains[[m]] <- advi
}



refaverage <- Reduce("+", listofrefchains)/nchains
newaverage <- Reduce("+", listofnewchains)/nchains
newstudentaverage <- Reduce("+", listofnewstudentchains)/nchains
malaaverage <- Reduce("+", listofmalachains)/nchains
adviaverage <- Reduce("+", listofadvichains)/nchains




#Autocorrelation
par(mfrow=c(1,6))
acf(refaverage[,1], main="RWM")
acf(newaverage[,1], main="IMH (Gaussian)")
acf(newstudentaverage[,1], main="IMH (Student)")
acf(malaaverage[,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(adviaverage[,1], main="ADVI")

par(mfrow=c(1,6))
acf(refaverage[,2], main="RWM")
acf(newaverage[,2], main="IMH (Gaussian)")
acf(newstudentaverage[,2], main="IMH (Student)")
acf(malaaverage[,2], main="MALA")
acf(vi[,2], main="NUTS")
acf(adviaverage[,2], main="ADVI")


par(mfrow=c(1,6))
acf(refaverage[,3], main="RWM")
acf(newaverage[,3], main="IMH (Gaussian)")
acf(newstudentaverage[,3], main="IMH (Student)")
acf(malaaverage[,3], main="MALA")
acf(vi[,3], main="NUTS")
acf(adviaverage[,3], main="ADVI")

#MSJD average
mssd(refaverage[,1])
mssd(newaverage[,1])
mssd(newstudentaverage[,1])
mssd(malaaverage[,1])
mssd(adviaverage[,1])
mssd(vi[,1])

mssd(refaverage[,2])
mssd(newaverage[,2])
mssd(newstudentaverage[,2])
mssd(malaaverage[,2])
mssd(adviaverage[,2])
mssd(vi[,2])

mssd(refaverage[,3])
mssd(newaverage[,3])
mssd(newstudentaverage[,3])
mssd(malaaverage[,3])
mssd(adviaverage[,3])
mssd(vi[,3])






#MSJD
mssd(ref[,1])
mssd(new[,1])
mssd(new.student[,1])
mssd(mala[,1])
mssd(advi[,1])
mssd(vi[,1])


#ESS
library("mcmcse")
ess(refaverage)
ess(newaverage)
ess(newstudentaverage)
ess(malaaverage)
ess(adviaverage)
ess(vi)


library(coda)
effectiveSize(refaverage)
effectiveSize(newaverage)
effectiveSize(newstudentaverage)
effectiveSize(malaaverage)
effectiveSize(adviaverage)
effectiveSize(vi)


effectiveSize(newaverage[1:20000,])

vi.obj <- as.mcmc(newaverage[1:2000,])
autocorr(vi.obj[,1]) 

vi.obj <- as.mcmc(newaverage[1:20000,])
autocorr(vi.obj[,1]) 



############# FOR ONE CHAIN ############



#Autocorrelation
par(mfrow=c(1,6))
acf(listofrefchains[[1]][,1], main="RWM")
acf(listofnewchains[[1]][,1], main="IMH (Gaussian)")
acf(listofnewstudentchains[[1]][,1], main="IMH (Student)")
acf(listofmalachains[[1]][,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(listofadvichains[[1]][,1], main="ADVI")

par(mfrow=c(1,6))
acf(listofrefchains[[1]][,2], main="RWM")
acf(listofnewchains[[1]][,2], main="IMH (Gaussian)")
acf(listofnewstudentchains[[1]][,2], main="IMH (Student)")
acf(listofmalachains[[1]][,2], main="MALA")
acf(vi[,2], main="NUTS")
acf(listofadvichains[[1]][,2], main="ADVI")


par(mfrow=c(1,6))
acf(listofrefchains[[1]][,3], main="RWM")
acf(listofnewchains[[1]][,3], main="IMH (Gaussian)")
acf(listofnewstudentchains[[1]][,3], main="IMH (Student)")
acf(listofmalachains[[1]][,3], main="MALA")
acf(vi[,3], main="NUTS")
acf(listofadvichains[[1]][,3], main="ADVI")

#MSJD
mssd(listofrefchains[[1]][,1])
mssd(listofnewchains[[1]][,1])
mssd(listofnewstudentchains[[1]][,1])
mssd(listofmalachains[[1]][,1])
mssd(listofadvichains[[1]][,1])
mssd(vi[,1])

mssd(listofrefchains[[1]][,2])
mssd(listofnewchains[[1]][,2])
mssd(listofnewstudentchains[[1]][,2])
mssd(listofmalachains[[1]][,2])
mssd(listofadvichains[[1]][,2])
mssd(vi[,2])

mssd(listofrefchains[[1]][,3])
mssd(listofnewchains[[1]][,3])
mssd(listofnewstudentchains[[1]][,3])
mssd(listofmalachains[[1]][,3])
mssd(listofadvichains[[1]][,3])
mssd(vi[,3])




#ESS
library("mcmcse")
ess(listofrefchains[[1]])
ess(listofnewchains[[1]])
ess(listofnewstudentchains[[1]])
ess(listofmalachains[[1]])
ess(listofadvichains[[1]])
ess(vi)


library(coda)
effectiveSize(listofrefchains[[1]])
effectiveSize(listofnewchains[[1]])
effectiveSize(listofnewstudentchains[[1]])
effectiveSize(listofmalachains[[1]])
effectiveSize(listofadvichains[[1]])
effectiveSize(vi)

