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


plotquantile1 <- function(df,df2,df3, title=NULL, ylim=NULL)
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


plotquantile2 <- function(df,df2,df3, title=NULL, ylim=NULL)
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


plotquantile3 <- function(df,df2,df3, title=NULL, ylim=NULL)
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
save1 = plotquantile1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantmala[,c(1,2,5)])
save2 = plotquantile2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantmala[,c(1,3,5)])
save3 = plotquantile3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantmala[,c(1,4,5)])

save1 = plotquantile1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantvi[,c(1,2,5)])
save2 = plotquantile2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantvi[,c(1,3,5)])
save3 = plotquantile3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantvi[,c(1,4,5)])

quantreftrue <- quantref
quantnewtrue <- quantnew
load("oldRData/quantvarinf.RData")

save1 = plotquantile1(quantreftrue[,c(1,2,5)],quantnewtrue[,c(1,2,5)],quantvarinf[,c(1,2,5)])
save2 = plotquantile2(quantreftrue[,c(1,3,5)],quantnewtrue[,c(1,3,5)],quantvarinf[,c(1,3,5)])
save3 = plotquantile3(quantreftrue[,c(1,4,5)],quantnewtrue[,c(1,4,5)],quantvarinf[,c(1,4,5)])

save <- grid.arrange(save1,save2,save3, ncol=3)

