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

load("rtte_mala.RData")



q1ref[,2] <- q1ref[,2] + 10
q2ref[,2] <- q2ref[,2] + 10
q3ref[,2] <- q3ref[,2] + 10

q1new[,2] <- q1new[,2] + 10
q2new[,2] <- q2new[,2] + 10
q3new[,2] <- q3new[,2] + 10


q1mala[,2] <- q1mala[,2] + 10
q2mala[,2] <- q2mala[,2] + 10
q3mala[,2] <- q3mala[,2] + 10



q1ref[,3] <- q1ref[,3] + 3 
q2ref[,3] <- q2ref[,3] + 3
q3ref[,3] <- q3ref[,3] + 3

q1new[,3] <- q1new[,3] + 3
q2new[,3] <- q2new[,3] + 3
q3new[,3] <- q3new[,3] + 3

q1mala[,3] <- q1mala[,3] + 3
q2mala[,3] <- q2mala[,3] + 3
q3mala[,3] <- q3mala[,3] + 3



quantref <- rbind(q1ref[-c(1:burn),],q2ref[-c(1:burn),],q3ref[-c(1:burn),])
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])
quantmala <- rbind(q1mala[-c(1:burn),],q2mala[-c(1:burn),],q3mala[-c(1:burn),])

colnames(quantref) <- colnames(quantnew)<-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")
colnames(quantmala) <-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")




plotquantile <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    if (j<3){
      grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(lambda[i])))
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    }else{
      grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(beta[i])))
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    }
    
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=2, top=title))
}


plotquant2<- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    if (j<3){
      grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(lambda[i])))
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    }else{
      grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=0.8) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 1,size=0.8)+
      xlab("iteration")+ theme_bw() +ylab(expression(paste(beta[i])))
      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2))
    }
    
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=2, top=title))
}


save = plotquantile(quantref,quantnew,quantmala)
save = grid.arrange(save,ncol=1)
ggsave(save, file="newpics/quantmala.pdf", width = 900, height = 450, units = "mm")

save = plotquant2(quantref,quantnew)
save = grid.arrange(save,ncol=1)
ggsave(save, file="newpics/quant_tte.pdf", width = 900, height = 450, units = "mm")


#########NUTS############
# q1nuts[,2] <- q1nuts[,2] + 10
# q2nuts[,2] <- q2nuts[,2] + 10
# q3nuts[,2] <- q3nuts[,2] + 10

# q1nuts[,3] <- q1nuts[,3] + 3
# q2nuts[,3] <- q2nuts[,3] + 3
# q3nuts[,3] <- q3nuts[,3] + 3

# quantnuts <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])
# colnames(quantnuts) <-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")

# save = plotquantile(quantref,quantnew,quantnuts)
# save = grid.arrange(save,ncol=1)
# ggsave(save, file="newpics/quantnuts.pdf", width = 900, height = 450, units = "mm")

#########ADVI############
# q1advi[,2] <- q1advi[,2] + 10
# q2advi[,2] <- q2advi[,2] + 10
# q3advi[,2] <- q3advi[,2] + 10

# q1advi[,3] <- q1advi[,3] + 3
# q2advi[,3] <- q2advi[,3] + 3
# q3advi[,3] <- q3advi[,3] + 3

# quantadvi <- rbind(q1advi[-c(1:burn),],q2advi[-c(1:burn),],q3advi[-c(1:burn),])
# colnames(quantadvi) <-c("iteration",expression(paste(lambda)),expression(paste(beta)),"quantile")

# save = plotquantile(quantref,quantnew,quantadvi)
# save = grid.arrange(save,ncol=1)
# ggsave(save, file="newpics/quantadvi.pdf", width = 900, height = 450, units = "mm")
