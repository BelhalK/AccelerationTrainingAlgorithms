library("mlxR")
library(saemix)


# ka_true <- 0.7272239
# V_true <- 7.536013
# k_true <- 1.756585
# o_ka <- sqrt(0.5078622)
# o_V <- sqrt(0.03980655)
# o_k <- sqrt(0.04624554)
# a_true<-1


# ka_true <- 1
# V_true <- 8
# k_true <- 2
# o_ka <- sqrt(0.1)
# o_V <- sqrt(0.05)
# o_k <- sqrt(0.07)
# a_true<-1


model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  Cl<-psi[id,3]
  k <- Cl/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

myModel <- inlineModel("


[INDIVIDUAL]
input = {ka_pop, V_pop, Cl_pop, omega_ka, omega_V, omega_Cl}
DEFINITION:
ka = {distribution=lognormal, reference=ka_pop, sd=omega_ka}
V  = {distribution=lognormal, reference=V_pop,  sd=omega_V }
Cl = {distribution=lognormal, reference=Cl_pop, sd=omega_Cl}


[LONGITUDINAL]
input = {ka, V, Cl,a}
EQUATION:
C = pkmodel(ka,V,Cl)
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=200

pop.param   <- c(
  ka_pop  = 1,    omega_ka  = 0.5,
  V_pop   = 10,   omega_V   = 0.4,
  Cl_pop  = 1,    omega_Cl  = 0.3, a =1)
  
res <- simulx(model     = myModel,
              parameter = pop.param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(0,24,by=0.5)))

 warfarin.saemix <- res$y
  warfarin.saemix$amount <- 100

  saemix.model<-saemixModel(model=model1cpt,description="warfarin"
  ,psi0=matrix(c(2,15,2,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","Cl"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")


K1 = 2000
K2 = 400
iterations = 1:(K1+K2+1)
end = K1+K2
options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),displayProgress=TRUE)
theo_ref<-saemix(saemix.model,saemix.data,options)



