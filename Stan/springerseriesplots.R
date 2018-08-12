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
# save.image("quantile_springer.RData")
load("quantile_springer.RData")

# quantref <- rbind(q1ref[-c(1:burn),],q2ref[-c(1:burn),],q3ref[-c(1:burn),])
# quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])
# quantmala <- rbind(q1mala[-c(1:burn),],q2mala[-c(1:burn),],q3mala[-c(1:burn),])
# quantnuts <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])

# quantadvi.full <- rbind(q1advi.full[-c(1:burn),],q2advi.full[-c(1:burn),],q3advi.full[-c(1:burn),])
# # quantadvi.onlymuvi <- rbind(q1advi.onlymuvi[-c(1:burn),],q2advi.onlymuvi[-c(1:burn),],q3advi.onlymuvi[-c(1:burn),])
# # quantadvi.onlygammavi <- rbind(q1advi.onlygammavi[-c(1:burn),],q2advi.onlygammavi[-c(1:burn),],q3advi.onlygammavi[-c(1:burn),])

# colnames(quantmala) <-colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
# colnames(quantnuts) <-c("iteration","ka","V","k","quantile")



plotq1 <- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",size=2, linetype=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",size=3, linetype="dotted")+
      xlab("iteration") + theme_bw() +ylab(expression(paste(ka[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=46), 
                 axis.title=element_text(size=48),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotq2 <- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",size=2, linetype=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",size=3, linetype="dotted")+
      xlab("iteration")+  theme_bw() +ylab(expression(paste(V[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=46), 
                 axis.title=element_text(size=48),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotq3 <- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",size=2, linetype=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",size=3, linetype="dotted")+
      xlab("iteration")+  theme_bw() +ylab(expression(paste(k[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=46), 
                 axis.title=element_text(size=48),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


iter = seq(1,59700, by=100)

colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
save1 = plotq1(quantref[iter,c(1,2,5)],quantnew[iter,c(1,2,5)])
save2 = plotq2(quantref[iter,c(1,3,5)],quantnew[iter,c(1,3,5)])
save3 = plotq3(quantref[iter,c(1,4,5)],quantnew[iter,c(1,4,5)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="baysmpics/quantref3.pdf", width = 900, height = 450, units = "mm")



plotq1_3 <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="red",size=3, linetype="dotted") +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="black",size=2, linetype="solid")+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="blue",size=2, linetype="dashed")+
      xlab("iteration")+  theme_bw() +ylab(expression(paste(ka[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=46), 
                 axis.title=element_text(size=48),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotq2_3 <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="red",size=3, linetype="dotted") +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="black",size=2, linetype="solid")+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="blue",size=2, linetype="dashed")+
      xlab("iteration")+  theme_bw() +ylab(expression(paste(V[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=46), 
                 axis.title=element_text(size=48),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotq3_3 <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="red",size=3, linetype="dotted") +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="black",size=2, linetype="solid")+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="blue",size=2, linetype="dashed")+
      xlab("iteration")+  theme_bw() +ylab(expression(paste(k[i])))
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=46), 
                 axis.title=element_text(size=48),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}




colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
save1 = plotq1_3(quantnew[iter,c(1,2,5)],quantmala[iter,c(1,2,5)],quantnuts[iter,c(1,2,5)])
save2 = plotq2_3(quantnew[iter,c(1,3,5)],quantmala[iter,c(1,3,5)],quantnuts[iter,c(1,3,5)])
save3 = plotq3_3(quantnew[iter,c(1,4,5)],quantmala[iter,c(1,4,5)],quantnuts[iter,c(1,4,5)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="baysmpics/quantmalaandnuts2.pdf", width = 900, height = 450, units = "mm")

