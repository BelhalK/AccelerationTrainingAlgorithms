library(saemix)
data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/annealing/data/joint.txt", header=T)
saemix.data<-saemixData(name.data=data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("dose","time"),name.response=c("y"), name.X="time")

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

saemix.model<-saemixModel(model=model1cpt,description="pk"
  ,psi0=matrix(c(1,10,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE), error.model="combined")

##RUNS

K1 = 300
K2 = 100
iterations = 1:(K1+K2)
end = K1+K2


#pkrin
options<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0)
pk<-saemix(saemix.model,saemix.data,options)


options.sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, alpha.sa=0.95)
pk.sa<-saemix(saemix.model,saemix.data,options.sa)
