library(saemix)


# zifro_data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/paramToRV/data/dataPK_zifrosilone.csv", header=T,sep=";")
zifro_data <- read.table("/Users/karimimohammedbelhal/Desktop/zifro/data/dataPK_zifrosilone.txt", header=T)
saemix.data_zifro<-saemixData(name.data=zifro_data,header=TRUE,sep=" ",na=NA, name.group=c("ID"),
  name.predictors=c("AMT","TIME"),name.response=c("Y"), name.X="X")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  T<-psi[id,1]
  ka<-psi[id,2]
  V<-psi[id,3]
  alpha<-psi[id,4]
  beta<-psi[id,5]
  CL<-alpha*V^beta
  k<-CL/V
  # ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*(tim-T))-exp(-ka*(tim-T)))
  return(ypred)
}


saemix.model_zifro<-saemixModel(model=model1cpt,description="zifrorin"
  ,psi0=matrix(c(0.2,1,250,1,1),ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
  transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5, 
  byrow=TRUE),error.model="exponential")

# saemix.model_zifronovar<-saemixModel(model=model1cpt,description="zifrorin"
#   ,psi0=matrix(c(0.2,1,250,1,1),ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
#   transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol=5, 
#   byrow=TRUE),error.model="exponential")

saemix.model_zifronovar<-saemixModel(model=model1cpt,description="zifrorin"
  ,psi0=matrix(c(0.158,0.18,40,1,1),ncol=5,byrow=TRUE, dimnames=list(NULL, c("T","ka","V","alpha","beta"))),
  transform.par=c(1,1,1,1,1),omega.init=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),ncol=5,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol=5, 
  byrow=TRUE),error.model="exponential")


K1 = 400
K2 = 200
iterations = 1:(K1+K2+1)
end = K1+K2


#With var no sa
options_zifro_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0)
zifro_without<-data.frame(saemix(saemix.model_zifro,saemix.data_zifro,options_zifro_without))
zifro_without <- cbind(iterations, zifro_without)

options_zifro_without<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2), nbiter.sa=0,nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0)
zifro_without<-data.frame(saemix(saemix.model_zifronovar,saemix.data_zifro,options_zifro_without))
zifro_without <- cbind(iterations, zifro_without)