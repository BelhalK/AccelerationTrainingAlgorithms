library("mlxR")
library(saemix)

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

model <- inlineModel("


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

N=5

param   <- c(
  ka_pop  = 1,    omega_ka  = 0.5,
  V_pop   = 10,   omega_V   = 0.4,
  Cl_pop  = 1,    omega_Cl  = 0.3, a =1)
  
res <- simulx(model     = model,
              parameter = param,
              treatment = list(time=0, amount=100),
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(1,3,by=1)))

 
res$y[2,3]
typeof(res$y[2,3])

writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Desktop/data_pk/res1a.txt", sep="\t")
b <- read.table("/Users/karimimohammedbelhal/Desktop/data_pk/res1a.txt", header=TRUE, sep="\t")
d <- readDatamlx(datafile = '/Users/karimimohammedbelhal/Desktop/data_pk/res1a.txt',  header   = c('id','time','y','amount'))

d$y[3,3]
typeof(d$y[2,3])

writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Desktop/data_pk/test.csv")
warfarin.saemix<-read.table("/Users/karimimohammedbelhal/Desktop/data_pk/test.csv", header=T, sep=",")
warfarin.saemix[3,3]
typeof(warfarin.saemix[2,3])


warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code/warfarin/warfarin_final.txt", header=T)
typeof(warfa_data[2,4])



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



