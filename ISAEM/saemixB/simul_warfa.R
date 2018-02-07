library("mlxR")
library(saemix)


ka_true <- 0.7272239
V_true <- 7.536013
k_true <- 1.756585
o_ka <- sqrt(0.5078622)
o_V <- sqrt(0.03980655)
o_k <- sqrt(0.04624554)
a_true<-1

model2 <- inlineModel("

                   [INDIVIDUAL]
                  input = {ka_pop, omega_ka, V_pop, omega_V, k_pop, omega_k}

                  DEFINITION:
                  ka = {distribution=lognormal, typical=ka_pop, sd=omega_ka}
                  V = {distribution=lognormal, typical=V_pop, sd=omega_V}
                  k = {distribution=lognormal, typical=k_pop, sd=omega_k} 
                  [LONGITUDINAL]
                  input = {ka,V, k, a}
                  EQUATION:
                  Cc = pkmodel(ka,V,k)
                  DEFINITION:
                  y = {distribution=normal, prediction=Cc, sd=a}
                          ")

    adm <- list(time=0, amount=100)
    y <- list(name='y', time=seq(1, 10, by=1))
    p <- c(ka_pop=ka_true, omega_ka=o_ka,
           V_pop=V_true, omega_V=o_V, 
           k_pop=k_true, omega_k=o_k,
           a=a_true)
    g <- list(size=100, level="individual")
    res<-simulx(model=model2, parameter=p, output=y, treatment=adm, group=g)
  
  warfarin.saemix <- res$y
  warfarin.saemix$amount <- 100

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}


saemix.model<-saemixModel(model=model1cpt,description="warfarin"
  ,psi0=matrix(c(3,10,2,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))


saemix.data<-saemixData(name.data=warfarin.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y"), name.X="time")


K1 = 600
K2 = 400
iterations = 1:(K1+K2+1)
end = K1+K2
options<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),displayProgress=TRUE)
theo_ref<-saemix(saemix.model,saemix.data,options)



