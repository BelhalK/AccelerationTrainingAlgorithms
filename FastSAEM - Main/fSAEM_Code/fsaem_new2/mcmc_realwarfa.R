load("newmcmc.RData")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
# save.image("newmcmc.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_new2/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('estep_mcmc.R')
  source('main.R')
  source('main_estep.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('mcmc.R')  
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  # source('SaemixRes_c.R') 
  source('SaemixObject.R') 
  source('zzz.R')
  
setwd("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_new2")
source('graphplot.R')


warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/research/CSDA/csda_new2/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

test = warfa_data[which(warfa_data[,2]>0),]
plotdata(test[,c(2,4,1)])
plotdata(warfa_data[,c(2,4,1)])

ggsave(file="pics/data.pdf", width = 900, height = 225, units = "mm")
ggsave(file="newpics/datap.pdf", width = 900, height = 450, units = "mm")

n <- length(unique(warfa_data$id))
model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]

  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

##RUNS

K1 = 400
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2

# #Warfarin
# options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
# warfa <- cbind(iterations, warfa)


# options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(K1))
# warfanew<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew))
# warfanew <- cbind(iterations, warfanew)


graphConvMC_twokernels(warfa,warfanew)


#compareMCMC

saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
  ,psi0=matrix(c(0.7,7.51,0.0178),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(0.2,0,0,0,0.18,0,0,0,0.03),ncol=3,byrow=TRUE),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE))

# saemix.model_warfa<-saemixModel(model=model1cpt,description="warfarin",type="structural"
#   ,psi0=matrix(c(0.5,7,0.1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
#   byrow=TRUE))


L_mcmc=10000
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta

options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,2,2,2),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta

# options_warfanew2<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,0,2,2,2),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# new2<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew2)$eta

# #Autocorrelation
# rwm.obj <- as.mcmc(ref[[10]])
# autocorr.plot(rwm.obj[,1]) + title("RWM SAEM Autocorrelation")

# new.obj <- as.mcmc(new[[10]])
# autocorr.plot(new.obj[,1]) + title("Laplace SAEM Autocorrelation")


# #MSJD
# mssd(ref[[10]][,1])
# mssd(new[[10]][,1])



#quantiles

i <- 10
qref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[[i]][1:k,dim], 0.2)
    qref[[dim]][k,2] <- quantile(ref[[i]][1:k,dim], 0.5)
    qref[[dim]][k,3] <- quantile(ref[[i]][1:k,dim], 0.8)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
for (dim in 1:3){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[[i]][1:k,dim], 0.2)
    qnew[[dim]][k,2] <- quantile(new[[i]][1:k,dim], 0.5)
    qnew[[dim]][k,3] <- quantile(new[[i]][1:k,dim], 0.8)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
}



plotquantile1 <- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(ka[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotquantile2 <- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(V[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotquantile3 <- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(k[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


iteration <- 1:L_mcmc
burn <- 100
q1ref <- data.frame(cbind(iteration,qref[[1]][,1],qref[[2]][,1],qref[[3]][,1]))
q2ref <- data.frame(cbind(iteration,qref[[1]][,2],qref[[2]][,2],qref[[3]][,2]))
q3ref <- data.frame(cbind(iteration,qref[[1]][,3],qref[[2]][,3],qref[[3]][,3]))
q1ref$quantile <- 1
q2ref$quantile <- 2
q3ref$quantile <- 3
quantref <- rbind(q1ref[-c(1:burn),],q2ref[-c(1:burn),],q3ref[-c(1:burn),])



q1new <- data.frame(cbind(iteration,qnew[[1]][,1],qnew[[2]][,1],qnew[[3]][,1]))
q2new <- data.frame(cbind(iteration,qnew[[1]][,2],qnew[[2]][,2],qnew[[3]][,2]))
q3new <- data.frame(cbind(iteration,qnew[[1]][,3],qnew[[2]][,3],qnew[[3]][,3]))
q1new$quantile <- 1
q2new$quantile <- 2
q3new$quantile <- 3
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])



q1ref[,2] <- q1ref[,2] + 1 
q2ref[,2] <- q2ref[,2] + 1
q3ref[,2] <- q3ref[,2] + 1

q1new[,2] <- q1new[,2] + 1
q2new[,2] <- q2new[,2] + 1
q3new[,2] <- q3new[,2] + 1


q1ref[,3] <- q1ref[,3] + 8 
q2ref[,3] <- q2ref[,3] + 8
q3ref[,3] <- q3ref[,3] + 8

q1new[,3] <- q1new[,3] + 8
q2new[,3] <- q2new[,3] + 8
q3new[,3] <- q3new[,3] + 8


q1ref[,4] <- q1ref[,4] + 0.01 
q2ref[,4] <- q2ref[,4] + 0.01
q3ref[,4] <- q3ref[,4] + 0.01

q1new[,4] <- q1new[,4] + 0.01
q2new[,4] <- q2new[,4] + 0.01
q3new[,4] <- q3new[,4] + 0.01


colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
save1 = plotquantile1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)])
save2 = plotquantile2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)])
save3 = plotquantile3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)])

save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save,file="pics_square/quantpknotlognew.pdf", width = 900, height = 225, units = "mm")

ggsave(save,file="newpics/quant_pk.pdf", width = 900, height = 450, units = "mm")

quantfinal <- rbind(quantref, quantnew)


dm1 <- melt(d[,1:n.var], id = NULL, value.name = "level", variable.name = "gene")
p1 <- ggplot(dm1, aes(x=gene, y=level)) + geom_boxplot(outlier.size=1.5) + 
  ylab("level")  
p1 <- p1 + theme(legend.position = "none", axis.text=element_text(size=18), 
                 axis.title=element_text(size=34),
                   axis.line = element_line(colour = 'black', size = 1.25),
                   axis.ticks = element_line(colour = "black", size = 1.25),
                   panel.border = element_blank())

i.figure <- i.figure + 1
png(figure.paper[i.figure], width = 900, height = 450)
print(p1)
dev.off()

