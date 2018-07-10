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
# load("hmc_quantile.RData")
load("hmc_quantile_indiv.RData")

l <- c(1,8,0.01)
for (d in 2:4){
  q1ref[,d] <- q1ref[,d] + l[d-1] 
  q2ref[,d] <- q2ref[,d] + l[d-1]
  q3ref[,d] <- q3ref[,d] + l[d-1]

  q1new[,d] <- q1new[,d] + l[d-1] 
  q2new[,d] <- q2new[,d] + l[d-1]
  q3new[,d] <- q3new[,d] + l[d-1]

  q1mala[,d] <- q1mala[,d] + l[d-1] 
  q2mala[,d] <- q2mala[,d] + l[d-1]
  q3mala[,d] <- q3mala[,d] + l[d-1]

  q1vi[,d] <- q1vi[,d] + l[d-1] 
  q2vi[,d] <- q2vi[,d] + l[d-1]
  q3vi[,d] <- q3vi[,d] + l[d-1]

  ###ADVI
  q1advi.full[,d] <- q1advi.full[,d] + l[d-1] 
  q2advi.full[,d] <- q2advi.full[,d] + l[d-1]
  q3advi.full[,d] <- q3advi.full[,d] + l[d-1]

  # q1advi.onlymuvi[,d] <- q1advi.onlymuvi[,d] + l[d-1] 
  # q2advi.onlymuvi[,d] <- q2advi.onlymuvi[,d] + l[d-1]
  # q3advi.onlymuvi[,d] <- q3advi.onlymuvi[,d] + l[d-1]

  # q1advi.onlygammavi[,d] <- q1advi.onlygammavi[,d] + l[d-1] 
  # q2advi.onlygammavi[,d] <- q2advi.onlygammavi[,d] + l[d-1]
  # q3advi.onlygammavi[,d] <- q3advi.onlygammavi[,d] + l[d-1]

}


quantref <- rbind(q1ref[-c(1:burn),],q2ref[-c(1:burn),],q3ref[-c(1:burn),])
quantnew <- rbind(q1new[-c(1:burn),],q2new[-c(1:burn),],q3new[-c(1:burn),])
quantmala <- rbind(q1mala[-c(1:burn),],q2mala[-c(1:burn),],q3mala[-c(1:burn),])
quantnuts <- rbind(q1vi[-c(1:burn),],q2vi[-c(1:burn),],q3vi[-c(1:burn),])

quantadvi.full <- rbind(q1advi.full[-c(1:burn),],q2advi.full[-c(1:burn),],q3advi.full[-c(1:burn),])
# quantadvi.onlymuvi <- rbind(q1advi.onlymuvi[-c(1:burn),],q2advi.onlymuvi[-c(1:burn),],q3advi.onlymuvi[-c(1:burn),])
# quantadvi.onlygammavi <- rbind(q1advi.onlygammavi[-c(1:burn),],q2advi.onlygammavi[-c(1:burn),],q3advi.onlygammavi[-c(1:burn),])

colnames(quantmala) <-colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
colnames(quantnuts) <-c("iteration","ka","V","k","quantile")

colnames(quantadvi.full) <-colnames(quantadvi.onlymuvi) <- colnames(quantadvi.onlygammavi)<-c("iteration","ka","V","k","quantile")
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
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
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
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
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
    grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}



colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantnew[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantnew[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantnew[,c(1,4,5)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="newpics/quant_pk.pdf", width = 900, height = 450, units = "mm")


colnames(quantref) <- colnames(quantnew)<-c("iteration","ka","V","k","quantile")
save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantmala[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantmala[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantmala[,c(1,4,5)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="newpics/quantmala.pdf", width = 900, height = 450, units = "mm")

save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantnuts[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantnuts[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantnuts[,c(1,4,5)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="newpics/quantnuts.pdf", width = 900, height = 450, units = "mm")


save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantadvi.full[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantadvi.full[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantadvi.full[,c(1,4,5)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="newpics/quantadvi_full.pdf", width = 900, height = 450, units = "mm")

save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantadvi.onlymuvi[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantadvi.onlymuvi[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantadvi.onlymuvi[,c(1,4,5)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="newpics/quantadvi_onlymuvi.pdf", width = 900, height = 450, units = "mm")

save1 = plotq1(quantref[,c(1,2,5)],quantnew[,c(1,2,5)],quantadvi.onlygammavi[,c(1,2,5)])
save2 = plotq2(quantref[,c(1,3,5)],quantnew[,c(1,3,5)],quantadvi.onlygammavi[,c(1,3,5)])
save3 = plotq3(quantref[,c(1,4,5)],quantnew[,c(1,4,5)],quantadvi.onlygammavi[,c(1,4,5)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="newpics/quantadvi_onlygammavi.pdf", width = 900, height = 450, units = "mm")



# quantreftrue <- quantref
# quantnewtrue <- quantnew
# load("oldRData/quantvarinf.RData")

# save1 = plotq1(quantreftrue[,c(1,2,5)],quantnewtrue[,c(1,2,5)],quantvarinf[,c(1,2,5)])
# save2 = plotq2(quantreftrue[,c(1,3,5)],quantnewtrue[,c(1,3,5)],quantvarinf[,c(1,3,5)])
# save3 = plotq3(quantreftrue[,c(1,4,5)],quantnewtrue[,c(1,4,5)],quantvarinf[,c(1,4,5)])

save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="newpics/quantadvi.pdf", width = 900, height = 450, units = "mm")
