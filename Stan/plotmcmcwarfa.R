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
load("hmc_quantile.RData")



q1ref[,2] <- q1ref[,2] + 1 
q2ref[,2] <- q2ref[,2] + 1
q3ref[,2] <- q3ref[,2] + 1

q1new[,2] <- q1new[,2] + 1
q2new[,2] <- q2new[,2] + 1
q3new[,2] <- q3new[,2] + 1

q1advi[,2] <- q1advi[,2] + 1
q2advi[,2] <- q2advi[,2] + 1
q3advi[,2] <- q3advi[,2] + 1

q1mala[,2] <- q1mala[,2] + 1
q2mala[,2] <- q2mala[,2] + 1
q3mala[,2] <- q3mala[,2] + 1

q1vi[,2] <- q1vi[,2] + 1
q2vi[,2] <- q2vi[,2] + 1
q3vi[,2] <- q3vi[,2] + 1


q1ref[,3] <- q1ref[,3] + 8 
q2ref[,3] <- q2ref[,3] + 8
q3ref[,3] <- q3ref[,3] + 8

q1new[,3] <- q1new[,3] + 8
q2new[,3] <- q2new[,3] + 8
q3new[,3] <- q3new[,3] + 8

q1advi[,3] <- q1advi[,3] + 8
q2advi[,3] <- q2advi[,3] + 8
q3advi[,3] <- q3advi[,3] + 8

q1mala[,3] <- q1mala[,3] + 8
q2mala[,3] <- q2mala[,3] + 8
q3mala[,3] <- q3mala[,3] + 8

q1vi[,3] <- q1vi[,3] + 8
q2vi[,3] <- q2vi[,3] + 8
q3vi[,3] <- q3vi[,3] + 8

q1ref[,4] <- q1ref[,4] + 0.01 
q2ref[,4] <- q2ref[,4] + 0.01
q3ref[,4] <- q3ref[,4] + 0.01

q1new[,4] <- q1new[,4] + 0.01
q2new[,4] <- q2new[,4] + 0.01
q3new[,4] <- q3new[,4] + 0.01

q1advi[,4] <- q1advi[,4] + 0.01
q2advi[,4] <- q2advi[,4] + 0.01
q3advi[,4] <- q3advi[,4] + 0.01

q1mala[,4] <- q1mala[,4] + 0.01
q2mala[,4] <- q2mala[,4] + 0.01
q3mala[,4] <- q3mala[,4] + 0.01

q1vi[,4] <- q1vi[,4] + 0.01
q2vi[,4] <- q2vi[,4] + 0.01
q3vi[,4] <- q3vi[,4] + 0.01


quantadvi <- rbind(q1advi[-c(1:burn),],q2advi[-c(1:burn),],q3advi[-c(1:burn),])
quantref <- rbind(q1ref[-c(1:burn),],q2ref[-c(1:burn),],q3ref[-c(1:burn),])
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])
quantmala <- rbind(q1mala[-c(1:burn),],q2mala[-c(1:burn),],q3mala[-c(1:burn),])
quantnuts <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])

colnames(quantadvi) <-colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
colnames(quantmala)<-colnames(quantnuts) <-c("iteration","ka","V","k","quantile")
plotq1 <- function(df,df2,df3, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  df3$quantile <- as.factor(df3$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(ka[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=24), 
                 axis.title=element_text(size=30),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotq2 <- function(df,df2,df3, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  df3$quantile <- as.factor(df3$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(V[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=24), 
                 axis.title=element_text(size=30),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotq3 <- function(df,df2,df3, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  df3$quantile <- as.factor(df3$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(k[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=24), 
                 axis.title=element_text(size=30),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}



colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantmala[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantmala[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantmala[,c(1,4,5)])

save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantvi[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantvi[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantvi[,c(1,4,5)])


save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantadvi[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantadvi[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantadvi[,c(1,4,5)])



quantreftrue <- quantref
quantnewtrue <- quantnew
load("oldRData/quantvarinf.RData")

save1 = plotq1(quantreftrue[,c(1,2,5)],quantnewtrue[,c(1,2,5)],quantvarinf[,c(1,2,5)])
save2 = plotq2(quantreftrue[,c(1,3,5)],quantnewtrue[,c(1,3,5)],quantvarinf[,c(1,3,5)])
save3 = plotq3(quantreftrue[,c(1,4,5)],quantnewtrue[,c(1,4,5)],quantvarinf[,c(1,4,5)])

save <- grid.arrange(save1,save2,save3, ncol=3)

