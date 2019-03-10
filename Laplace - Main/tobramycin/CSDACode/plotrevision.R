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
 
 d <- ncol(ref)
i <- 10
start_interval <- 200
zero <- as.data.frame(matrix(0,nrow = L_mcmc-start_interval,ncol = 3))


#quantiles
qlow <- 0.2
qmed <- 0.5
qhigh <- 0.8


qref <- list(ref[1:L_mcmc,],ref[1:L_mcmc,],ref[1:L_mcmc,])
for (dim in 1:d){
  print(dim)
  for (k in 1:L_mcmc){
    qref[[dim]][k,1] <- quantile(ref[1:k,dim], qlow)
    qref[[dim]][k,2] <- quantile(ref[1:k,dim], qmed)
    qref[[dim]][k,3] <- quantile(ref[1:k,dim], qhigh)
  }
  qref[[dim]]$iteration <- 1:L_mcmc
}


qnew <- list(new[1:L_mcmc,],new[1:L_mcmc,],new[1:L_mcmc,])
for (dim in 1:d){
  print(dim)
  for (k in 1:L_mcmc){
    qnew[[dim]][k,1] <- quantile(new[1:k,dim], qlow)
    qnew[[dim]][k,2] <- quantile(new[1:k,dim], qmed)
    qnew[[dim]][k,3] <- quantile(new[1:k,dim], qhigh)
  }
  qnew[[dim]]$iteration <- 1:L_mcmc
}


qnew.student <- list(new.student[1:L_mcmc,],new.student[1:L_mcmc,],new.student[1:L_mcmc,])
for (dim in 1:d){
  print(dim)
  for (k in 1:L_mcmc){
    qnew.student[[dim]][k,1] <- quantile(new.student[1:k,dim], qlow)
    qnew.student[[dim]][k,2] <- quantile(new.student[1:k,dim], qmed)
    qnew.student[[dim]][k,3] <- quantile(new.student[1:k,dim], qhigh)
  }
  qnew.student[[dim]]$iteration <- 1:L_mcmc
}



iteration <- 1:L_mcmc
burn <- 100

q1ref <- data.frame(cbind(iteration,qref[[1]][,1],qref[[2]][,1]))
q2ref <- data.frame(cbind(iteration,qref[[1]][,2],qref[[2]][,2]))
q3ref <- data.frame(cbind(iteration,qref[[1]][,3],qref[[2]][,3]))
q1ref$quantile <- 1
q2ref$quantile <- 2
q3ref$quantile <- 3
quantref <- rbind(q2ref[-c(1:burn),])
colnames(quantref) <- c("iteration","V","k","quantile")

q1new <- data.frame(cbind(iteration,qnew[[1]][,1],qnew[[2]][,1]))
q2new <- data.frame(cbind(iteration,qnew[[1]][,2],qnew[[2]][,2]))
q3new <- data.frame(cbind(iteration,qnew[[1]][,3],qnew[[2]][,3]))
q1new$quantile <- 1
q2new$quantile <- 2
q3new$quantile <- 3
quantnew <- rbind(q2new[-c(1:burn),])
colnames(quantnew)<-c("iteration","V","k","quantile")


q1new.student <- data.frame(cbind(iteration,qnew.student[[1]][,1],qnew.student[[2]][,1]))
q2new.student <- data.frame(cbind(iteration,qnew.student[[1]][,2],qnew.student[[2]][,2]))
q3new.student <- data.frame(cbind(iteration,qnew.student[[1]][,3],qnew.student[[2]][,3]))
q1new.student$quantile <- 1
q2new.student$quantile <- 2
q3new.student$quantile <- 3
quantnew.student <- rbind(q2new.student[-c(1:burn),])
colnames(quantnew.student)<-c("iteration","V","k","quantile")

plotquantile(quantref,quantnew)


plotquantile.univariate1 <- function(df,df2, title=NULL, ylim=NULL)
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
      xlab("")+ylab(expression(paste(ka[i])))

      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
     if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotquantile.univariate2 <- function(df,df2, title=NULL, ylim=NULL)
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
      xlab("")+ylab(expression(paste(V[i])))

      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
     if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


plotquantile.univariate3 <- function(df,df2, title=NULL, ylim=NULL)
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
      xlab("")+ylab(expression(paste(k[i])))

      grafj <- grafj + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=34), 
                 axis.title=element_text(size=40),
                   panel.border = element_rect(colour = "black", fill=NA, size=2),plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"))
     if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}



## MEDIAN REF VS LAPLACE
save1 = plotquantile.univariate1(quantref[,c(1,2,4)],quantnew[,c(1,2,4)])
save2 = plotquantile.univariate2(quantref[,c(1,3,4)],quantnew[,c(1,3,4)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save, file="revisionpics/median_pk.pdf", width = 900, height = 300, units = "mm")

ggsave(save1, file="pics_newlayout/median_pk_ka.pdf", width = 900, height = 250, units = "mm")
ggsave(save2, file="pics_newlayout/median_pk_V.pdf", width = 900, height = 250, units = "mm")
ggsave(save3, file="pics_newlayout/median_pk_k.pdf", width = 900, height = 250, units = "mm")





## MEDIAN REF VS LAPLACE
save1 = plotquantile.univariate1.df3(quantref[,c(1,2,4)],quantnew[,c(1,2,4)],quantnew.student[,c(1,2,4)])
save2 = plotquantile.univariate2.df3(quantref[,c(1,3,4)],quantnew[,c(1,3,4)],quantnew.student[,c(1,3,4)])
save <- grid.arrange(save1,save2,save3, ncol=3)
ggsave(save1, file="pics_newlayout/median_pk_ka_student.pdf", width = 900, height = 250, units = "mm")
ggsave(save2, file="pics_newlayout/median_pk_V_student.pdf", width = 900, height = 250, units = "mm")
ggsave(save3, file="pics_newlayout/median_pk_k_student.pdf", width = 900, height = 250, units = "mm")



par(mfrow=c(2,3))
acf(ref[,1], main="RWM")
acf(ref[,2], main="RWM")
acf(ref[,3], main="RWM")
acf(new[,1], main="IMH")
acf(new[,2], main="IMH")
acf(new[,3], main="IMH")


png(filename="revisionpics/acf_two.png", units="in", width=10, height=10, res=300)
par(mfrow=c(3,2))
acf(newaverage[,1], main="IMH")
acf(malaaverage[,1], main="RWM")
acf(newaverage[,2], main="")
acf(malaaverage[,2], main="")
acf(newaverage[,3], main="")
acf(malaaverage[,3], main="")

dev.off()

png(filename="revisionpics/acf_two_student.png", units="in", width=10, height=3, res=300)
par(mfrow=c(1,3))
acf(newaverage[,1], main="IMH (Gaussian)")
acf(newstudentaverage[,1], main="IMH (Student)")
acf(malaaverage[,1], main="RWM")
dev.off()


png(filename="pics_newlayout/acf_two_ka.png", units="in", width=10, height=4, res=300)
par(mfrow=c(1,2))
acf(newaverage[,1], main="IMH")
acf(malaaverage[,1], main="RWM")
dev.off()

png(filename="pics_newlayout/acf_two_V.png", units="in", width=10, height=4, res=300)
par(mfrow=c(1,2))
acf(newaverage[,2], main="IMH")
acf(malaaverage[,2], main="RWM")
dev.off()

png(filename="pics_newlayout/acf_two_k.png", units="in", width=10, height=4, res=300)
par(mfrow=c(1,2))
acf(newaverage[,3], main="IMH")
acf(malaaverage[,3], main="RWM")
dev.off()




#####STUDENT


png(filename="pics_newlayout/acf_two_ka_student.png", units="in", width=10, height=4, res=300)
par(mfrow=c(1,3))
acf(newaverage[,1], main="IMH (Gaussian)")
acf(newstudentaverage[,1], main="IMH (Student)")
acf(malaaverage[,1], main="RWM")
dev.off()

png(filename="pics_newlayout/acf_two_V_student.png", units="in", width=10, height=4, res=300)
par(mfrow=c(1,3))
acf(newaverage[,2], main="IMH (Gaussian)")
acf(newstudentaverage[,2], main="IMH (Student)")
acf(malaaverage[,2], main="RWM")
dev.off()

png(filename="pics_newlayout/acf_two_k_student.png", units="in", width=10, height=4, res=300)
par(mfrow=c(1,3))
acf(newaverage[,3], main="IMH (Gaussian)")
acf(newstudentaverage[,3], main="IMH (Student)")
acf(malaaverage[,3], main="RWM")
dev.off()

### ACF ALL METHODS


#Autocorrelation
png(filename="revisionpics/acf_all.png", units="in", width=10, height=3, res=300)
par(mfrow=c(1,6))
acf(newaverage[,1], main="IMH (Gaussian)")
acf(newstudentaverage[,1], main="IMH (Student)")
acf(malaaverage[,1], main="RWM")
acf(refaverage[,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(adviaverage[,1], main="ADVI")
dev.off()




#Autocorrelation
png(filename="pics_newlayout/acf_all.png", units="in", width=12, height=2.6, res=300)
par(mfrow=c(1,5))
acf(newaverage[,1], main="IMH")
acf(malaaverage[,1], main="RWM")
acf(refaverage[,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(adviaverage[,1], main="ADVI")
dev.off()




#Autocorrelation
png(filename="pics_newlayout/acf_all_student.png", units="in", width=12, height=2.6, res=300)
par(mfrow=c(1,6))
acf(newaverage[,1], main="IMH (Gaussian)")
acf(newstudentaverage[,1], main="IMH (Student)")
acf(malaaverage[,1], main="RWM")
acf(refaverage[,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(adviaverage[,1], main="ADVI")
dev.off()



#Autocorrelation
png(filename="revisionpics/acf_all_layout.png", units="in", width=10, height=8, res=300)
par(mfrow=c(2,3))
acf(newaverage[,1], main="IMH (Gaussian)")
acf(newstudentaverage[,1], main="IMH (Student)")
acf(malaaverage[,1], main="RWM")
acf(refaverage[,1], main="MALA")
acf(vi[,1], main="NUTS")
acf(adviaverage[,1], main="ADVI")
dev.off()

# par(mfrow=c(1,6))
# acf(refaverage[,2], main="RWM")
# acf(newaverage[,2], main="IMH (Gaussian)")
# acf(newstudentaverage[,2], main="IMH (Student)")
# acf(malaaverage[,2], main="MALA")
# acf(vi[,2], main="NUTS")
# acf(adviaverage[,2], main="ADVI")


# par(mfrow=c(1,6))
# acf(refaverage[,3], main="RWM")
# acf(newaverage[,3], main="IMH (Gaussian)")
# acf(newstudentaverage[,3], main="IMH (Student)")
# acf(malaaverage[,3], main="MALA")
# acf(vi[,3], main="NUTS")
# acf(adviaverage[,3], main="ADVI")



## TABLES


#MSJD average
mssd(refaverage)
mssd(newaverage)
mssd(newstudentaverage)
mssd(malaaverage)
mssd(adviaverage)
mssd(vi)


#MSJD
mssd(ref)
mssd(new)
mssd(new.student)
mssd(mala)
mssd(advi)
mssd(vi)


#ESS
library("mcmcse")
ess(ref)
ess(new)
ess(new.student)
ess(mala)
ess(advi)
ess(vi)

#ESS
library("mcmcse")
ess(refaverage)
ess(newaverage)
ess(newstudentaverage)
ess(malaaverage)
ess(adviaverage)
ess(vi)


