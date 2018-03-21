require(ggplot2)
require(gridExtra)
require(reshape2)
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/annealing/R")
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

setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/annealing")
source('plots.R') 
# data <- read.csv("/Users/karimimohammedbelhal/Documents/GitHub/saem/annealing/data/joint_datanew.csv", header=T,sep=",")
# data <- data[data$ytype==1,]

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

saemix.model<-saemixModel(model=model1cpt,description="pk",type="structural"
  ,psi0=matrix(c(1,10,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE), error.model="combined")

##RUNS

K1 = 1000
K2 = 100
iterations = 1:(K1+K2)
end = K1+K2


#pkrin
options<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1)
pk<-data.frame(saemix(saemix.model,saemix.data,options))
pk <- cbind(iterations, pk[-1,])

options.sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),displayProgress=TRUE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1, alpha.sa=0.95)
pk.sa<-data.frame(saemix(saemix.model,saemix.data,options.sa))
pk.sa <- cbind(iterations, pk.sa[-1,])

plot2(pk,pk.sa)

optionsnew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, an=TRUE,coeff=2)
pknew<-data.frame(saemix(saemix.model,saemix.data,optionsnew))
pknew <- cbind(iterations, pknew[-1,])


plot3(pk,pk.sa,pknew)

optionsnew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, an=TRUE,coeff=2.5)
pknew2<-data.frame(saemix(saemix.model,saemix.data,optionsnew))
pknew2 <- cbind(iterations, pknew2[-1,])

optionsnew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, an=TRUE,coeff=1.5)
pknew3<-data.frame(saemix(saemix.model,saemix.data,optionsnew))
pknew3 <- cbind(iterations, pknew3[-1,])
plot5(pk,pk.sa,pknew,pknew2,pknew3)


replicate = 3

final.ref <- 0
final.sa <- 0
final.an <- 0
for (m in 1:replicate){
  print(m)
  l = list(c(1,10,1,0,0,0),c(2,5,2,0,0,0),c(3,12,3,0,0,0))
  
  saemix.model<-saemixModel(model=model1cpt,description="pk",type="structural"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE), error.model="combined")

  options<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1)
  pk<-data.frame(saemix(saemix.model,saemix.data,options))
  pk <- cbind(iterations, pk[-1,])
  pk['individual'] <- m
  final.ref <- rbind(final.ref,pk)


  options.sa<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),displayProgress=FALSE,nbiter.burn =0, map.range=c(0),an=FALSE,coeff=1, alpha.sa=0.95)
  pk.sa<-data.frame(saemix(saemix.model,saemix.data,options.sa))
  pk.sa <- cbind(iterations, pk.sa[-1,])
  pk.sa['individual'] <- m
  final.sa <- rbind(final.sa,pk.sa)

  optionsnew<-list(seed=39546,map=F,fim=F,ll.is=T,nbiter.mcmc = c(2,2,2,0), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE,nbiter.burn =0, an=TRUE,coeff=2)
  pknew<-data.frame(saemix(saemix.model,saemix.data,optionsnew))
  pknew <- cbind(iterations, pknew[-1,])
  pknew['individual'] <- m
  final.an <- rbind(final.an,pknew)

}




diff <- function(df,df2,df3, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  df3$individual <- as.factor(df3$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",size=1) +
    geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue",linetype = 2,size=1)+
    geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red",linetype = 2,size=1)+
      xlab("") +scale_x_log10(breaks= c(10,100,200))+ ylab(names(df[j])) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text( color="black", 
                           size=10, angle=0),
          axis.text.y = element_text( color="black", 
                           size=10, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black",  size=10)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}

diff(final.ref,final.sa,final.an)


