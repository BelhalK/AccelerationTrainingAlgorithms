load("realwarfa_mcmc_conv.RData")
# load("realwarfa_mcmc_newquantiles.RData")
# load("realwarfa_mcmc_conv_new.RData")
library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
# load("realwarfa_mcmc_conv_var.RData")
# save.image("realwarfa_mcmc_conv_varwithnew.RData")
# save.image("realwarfa_mcmc_newquantiles.RData")
# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
setwd("/Users/karimimohammedbelhal/Desktop/ongoing_research/CSDA/csda_new/R")
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
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  # source('SaemixRes_c.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  
setwd("/Users/karimimohammedbelhal/Desktop/ongoing_research/CSDA/csda_new")
source('graphplot.R')


warfa_data <- read.table("/Users/karimimohammedbelhal/Desktop/ongoing_research/CSDA/csda_new/data/warfarin_data.txt", header=T)
saemix.data_warfa<-saemixData(name.data=warfa_data,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")


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

#Warfarin
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
warfa<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfa))
warfa <- cbind(iterations, warfa)


options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,nbiter.mcmc = c(2,2,2,6), nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(K1))
warfanew<-data.frame(saemix(saemix.model_warfa,saemix.data_warfa,options_warfanew))
warfanew <- cbind(iterations, warfanew)


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


L_mcmc=3000
options_warfa<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
ref<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta_ref

options_warfanew<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
new<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta

# options_warfanew2<-list(seed=39546,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,0,0,2,2,2),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# new2<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew2)$eta

# #Autocorrelation
# rwm.obj <- as.mcmc(ref[[10]])
# autocorr.plot(rwm.obj[,1]) + title("RWM SAEM Autocorrelation")

# new.obj <- as.mcmc(new[[10]])
# autocorr.plot(new.obj[,1]) + title("Laplace SAEM Autocorrelation")

# new2.obj <- as.mcmc(new2[[10]])
# autocorr.plot(new2.obj[,1]) + title("new2 Autocorrelation")

# #MSJD
# mssd(ref[[10]][,1])
# mssd(new[[10]][,1])
# mssd(new2[[10]][,1])


# start_interval <- 200
# zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 3))

# etabarref <- 1/n*Reduce("+",ref)
# expecref <- data.frame(apply(etabarref[-(1:start_interval),], 2, cummean))
# expecref$iteration <- 1:(L_mcmc-start_interval)


# sdref <- 0
# for (i in 1:n){
#   var <- data.frame(apply(ref[[i]][-(1:start_interval),]^2, 2, cummean))
#   meansq <- data.frame(apply(ref[[i]][-(1:start_interval),], 2, cummean))^2
#   sdref <- sdref + sqrt(pmax(zero,var - meansq))
# }


# sdref <- 1/n*sdref
# sdref$iteration <- 1:(L_mcmc-start_interval)


# etabarnew <- 1/n*Reduce("+",new)
# expecnew <- data.frame(apply(etabarnew[-(1:start_interval),], 2, cummean))
# expecnew$iteration <- 1:(L_mcmc-start_interval)


# sdnew <- 0
# for (i in 1:n){
#   var <- data.frame(apply(new[[i]][-(1:start_interval),]^2, 2, cummean))
#   meansq <- data.frame(apply(new[[i]][-(1:start_interval),], 2, cummean))^2
#   sdnew <- sdnew + sqrt(pmax(zero,var - meansq))
# }

# sdnew <- 1/n*sdnew
# sdnew$iteration <- 1:(L_mcmc-start_interval)


# plotmcmc(expecref[,c(4,1:3)],expecnew[,c(4,1:3)],title="mean")
# plotmcmc(sdref[-c(1:10),c(4,1:3)],sdnew[-c(1:10),c(4,1:3)],title="sd")



# #one invdiv


# start_interval <- 200
# zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 3))


# for (i in 1:2){
# indetabarref <- ref[[i]]
# indexpecref <- data.frame(apply(indetabarref[-(1:start_interval),], 2, cummean))
# indexpecref$iteration <- 1:(L_mcmc-start_interval)


# indsdref <- 0
# indvar <- data.frame(apply(ref[[i]][-(1:start_interval),]^2, 2, cummean))
# indmeansq <- data.frame(apply(ref[[i]][-(1:start_interval),], 2, cummean))^2
# indsdref <- indsdref + sqrt(pmax(zero,indvar - indmeansq))
# indsdref$iteration <- 1:(L_mcmc-start_interval)


# indetabarnew <- new[[i]]
# indexpecnew <- data.frame(apply(indetabarnew[-(1:start_interval),], 2, cummean))
# indexpecnew$iteration <- 1:(L_mcmc-start_interval)


# indsdnew <- 0
# indvar <- data.frame(apply(new[[i]][-(1:start_interval),]^2, 2, cummean))
# indmeansq <- data.frame(apply(new[[i]][-(1:start_interval),], 2, cummean))^2
# indsdnew <- indsdnew + sqrt(pmax(zero,indvar - indmeansq))
# indsdnew$iteration <- 1:(L_mcmc-start_interval)

# indetabarnew2 <- new2[[i]]
# indexpecnew2 <- data.frame(apply(indetabarnew2[-(1:start_interval),], 2, cummean))
# indexpecnew2$iteration <- 1:(L_mcmc-start_interval)


# indsdnew2 <- 0
# indvar <- data.frame(apply(new2[[i]][-(1:start_interval),]^2, 2, cummean))
# indmeansq <- data.frame(apply(new2[[i]][-(1:start_interval),], 2, cummean))^2
# indsdnew2 <- indsdnew2 + sqrt(pmax(zero,indvar - indmeansq))
# indsdnew2$iteration <- 1:(L_mcmc-start_interval)

# plotmcmc(indexpecref[,c(4,1:3)],indexpecnew[,c(4,1:3)],title=paste("mean",i))
# plotmcmc(indsdref[-c(1:10),c(4,1:3)],indsdnew[-c(1:10),c(4,1:3)],title=paste("sd",i))
# }


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
  # plotmcmc(qref[[dim]][,c(4,1:3)],qnew[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
}

# for (dim in 1:3){
# plotmcmc(qref[[dim]][,c(4,1:3)],qnew[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
# }


# qnew2 <- list(new2[[i]][1:L_mcmc,],new2[[i]][1:L_mcmc,],new2[[i]][1:L_mcmc,])
# for (dim in 1:3){
#   print(dim)
#   for (k in 1:L_mcmc){
#     qnew2[[dim]][k,1] <- quantile(new2[[i]][1:k,dim], 0.05)
#     qnew2[[dim]][k,2] <- quantile(new2[[i]][1:k,dim], 0.5)
#     qnew2[[dim]][k,3] <- quantile(new2[[i]][1:k,dim], 0.95)
#   }
#   qnew2[[dim]]$iteration <- 1:L_mcmc
#   # plotmcmc(qref[[dim]][,c(4,1:3)],qnew2[[dim]][,c(4,1:3)],title=paste("quantiles",i,"dim", dim))
# }



plotquantile <- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=1)+
      xlab("")+scale_x_log10()+ theme_bw() +ylab(names(df[j]))+ theme(axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=15, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=20)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
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



colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
plotquantile(quantref,quantnew)

quantfinal <- rbind(quantref, quantnew)


dm1 <- melt(d[,1:n.var], id = NULL, value.name = "level", variable.name = "gene")
p1 <- ggplot(dm1, aes(x=gene, y=level)) + geom_boxplot(outlier.size=1.5) + 
  ylab("level")  
p1 <- p1 + theme(legend.position = "none", axis.text=element_text(size=18), 
                 axis.title=element_text(size=24),
                   axis.line = element_line(colour = 'black', size = 1.25),
                   axis.ticks = element_line(colour = "black", size = 1.25),
                   panel.border = element_blank())

i.figure <- i.figure + 1
png(figure.paper[i.figure], width = 900, height = 450)
print(p1)
dev.off()



# q1new2 <- data.frame(cbind(iteration,qnew2[[1]][,1],qnew2[[2]][,1],qnew2[[3]][,1]))
# q2new2 <- data.frame(cbind(iteration,qnew2[[1]][,2],qnew2[[2]][,2],qnew2[[3]][,2]))
# q3new2 <- data.frame(cbind(iteration,qnew2[[1]][,3],qnew2[[2]][,3],qnew2[[3]][,3]))
# q1new2$quantile <- 1
# q2new2$quantile <- 2
# q3new2$quantile <- 3
# quantnew2 <- rbind(q1new2[-c(1:burn),],q2new2[-c(1:burn),],q3new2[-c(1:burn),])


# colnames(quantnew2)<-c("iteration","ka","V","k","quantile")

# plotquantile3 <- function(df,df2,df3, title=NULL, ylim=NULL)
# {
#  G <- (ncol(df)-2)/3
#   df$quantile <- as.factor(df$quantile)
#   df2$quantile <- as.factor(df2$quantile)
#   df3$quantile <- as.factor(df3$quantile)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=1)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype = 2,size=1)+
#       xlab("")+scale_x_log10()+ theme_bw() +ylab(names(df[j]))+ theme(axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
#                            size=15, angle=0),
#           axis.text.y = element_text(face="bold", color="black", 
#                            size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=20)) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj

#   }
#   do.call("grid.arrange", c(graf, ncol=3, top=title))
# }

# plotquantile3(quantref[1:29700,],quantnew[1:29700,],quantnew2)
# # gelman.plot(mcmc.list(as.mcmc(ref[[10]])), bin.width = 10, max.bins = 50,confidence = 0.95, transform = FALSE, autoburnin=TRUE, auto.layout = TRUE)

# geweke.plot(mcmc.list(as.mcmc(ref[[10]])), frac1=0.1, frac2=0.5)
# geweke.plot(mcmc.list(as.mcmc(new[[10]])), frac1=0.1, frac2=0.5)
# #cdf

# i <- 10
# cdfref <- list(ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,],ref[[i]][1:L_mcmc,])

# for (dim in 1:3){
#   print(dim)
#   qf1 <- quantile(ref[[i]][,dim], 0.05)
#   qf2 <- quantile(ref[[i]][,dim], 0.5)
#   qf3 <- quantile(ref[[i]][,dim], 0.95)

#   for (k in 1:L_mcmc){
#     cdfref[[dim]][k,1] <- mean(ref[[i]][which(ref[[i]][1:k,dim] < qf1),dim])
#     cdfref[[dim]][k,2] <- mean(ref[[i]][which(ref[[i]][1:k,dim] < qf2),dim])
#     cdfref[[dim]][k,3] <- mean(ref[[i]][which(ref[[i]][1:k,dim] < qf3),dim])
#   }
#   cdfref[[dim]]$iteration <- 1:L_mcmc
# }


# cdfnew <- list(new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,],new[[i]][1:L_mcmc,])
# for (dim in 1:3){
#   print(dim)
#   qf1 <- quantile(new[[i]][,dim], 0.05)
#   qf2 <- quantile(new[[i]][,dim], 0.5)
#   qf3 <- quantile(new[[i]][,dim], 0.95)

#   for (k in 1:L_mcmc){
#     cdfnew[[dim]][k,1] <- mean(new[[i]][which(new[[i]][1:k,dim] < qf1),dim])
#     cdfnew[[dim]][k,2] <- mean(new[[i]][which(new[[i]][1:k,dim] < qf2),dim])
#     cdfnew[[dim]][k,3] <- mean(new[[i]][which(new[[i]][1:k,dim] < qf3),dim])
#   }
#   cdfnew[[dim]]$iteration <- 1:L_mcmc
# }



# iteration <- 1:L_mcmc

# cdf1ref <- data.frame(cbind(iteration,cdfref[[1]][,1],cdfref[[2]][,1],cdfref[[3]][,1]))
# cdf2ref <- data.frame(cbind(iteration,cdfref[[1]][,2],cdfref[[2]][,2],cdfref[[3]][,2]))
# cdf3ref <- data.frame(cbind(iteration,cdfref[[1]][,3],cdfref[[2]][,3],cdfref[[3]][,3]))
# cdf1ref$quantile <- 1
# cdf2ref$quantile <- 2
# cdf3ref$quantile <- 3
# cdfref <- rbind(cdf1ref[-c(1:10),],cdf2ref[-c(1:10),],cdf3ref[-c(1:10),])



# cdf1new <- data.frame(cbind(iteration,cdfnew[[1]][,1],cdfnew[[2]][,1],cdfnew[[3]][,1]))
# cdf2new <- data.frame(cbind(iteration,cdfnew[[1]][,2],cdfnew[[2]][,2],cdfnew[[3]][,2]))
# cdf3new <- data.frame(cbind(iteration,cdfnew[[1]][,3],cdfnew[[2]][,3],cdfnew[[3]][,3]))
# cdf1new$quantile <- 1
# cdf2new$quantile <- 2
# cdf3new$quantile <- 3
# cdfnew <- rbind(cdf1new[-c(1:10),],cdf2new[-c(1:10),],cdf3new[-c(1:10),])


# plotquantile(cdfref,cdfnew, title= "cdf")


# plotmcmc(indexpecref[,c(4,1:3)],indexpecnew[,c(4,1:3)],title=paste("mean",i))
# plotmcmc(indsdref[-c(1:10),c(4,1:3)],indsdnew[-c(1:10),c(4,1:3)],title=paste("sd",i))


# #asymptotic variance
# seed0 <- 39546
# replicate <- 5
# L_mcmc<-5000
# final.ref <- 0
# final.new <- 0
# for (m in  1:replicate){
# options_warfa<-list(seed=m*seed0,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(2,2,2,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# rwm<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfa)$eta_ref[[10]]
# rwm['individual'] <- m
# final.ref <- rbind(final.ref,rwm)

# options_warfanew<-list(seed=m*seed0,map=F,fim=F,ll.is=F,L_mcmc=L_mcmc,nbiter.mcmc = c(0,0,0,6,0),nb.chains=1, nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE,nbiter.burn =0, map.range=c(0))
# laplace<-mcmc(saemix.model_warfa,saemix.data_warfa,options_warfanew)$eta[[10]]
# laplace['individual'] <- m
# final.new <- rbind(final.new,laplace)
# }


# sd.ref <- 0
# sd.new <- 0
# for (m in  1:replicate){
#   mean <- data.frame(apply(final.ref[final.ref$individual==m,1:3], 2, cummean))
#   ML <- final.ref[final.ref$individual==m,1:3]
#   ML[1:(L_mcmc),]<- mean[L_mcmc,1:3]
#   sd.ref <- sd.ref+(mean - ML)^2

#   mean <- data.frame(apply(final.new[final.new$individual==m,1:3], 2, cummean))
#   ML <- final.new[final.new$individual==m,1:3]
#   ML[1:(L_mcmc),]<- mean[L_mcmc,1:3]
#   sd.new <- sd.new+(mean - ML)^2

# }

# sd.ref <- 1/replicate*sd.ref
# sd.new <- 1/replicate*sd.new
# sd.ref <- data.frame(sd.ref)
# sd.new <- data.frame(sd.new)

# sd.ref$iteration <- 1:L_mcmc
# sd.new$iteration <- 1:L_mcmc


# plotmcmc(sd.ref[1:2000,c(4,1:3)],sd.new[1:2000,c(4,1:3)],title="as var")



