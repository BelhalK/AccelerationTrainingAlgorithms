# require(ggplot2)
# require(gridExtra)
# require(reshape2)
library(rlist)
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
library(data.table)
source("utils/algos.R")
source("utils/func.R")
options(digits = 22)
# load("precisionagainstn_smalln.RData")


eml <- ieml <- iemseql <- oemvrl <- sagal <- list()
emiter <- iemiter <-iemseqiter <- oemvriter <- sagaiter <- list()

# datasizes <- c(seq(1000, 9000, 1000),seq(10000, 106000, 5000))
datasizes <- seq(1000, 10000, 1000)
nsim=5

emindex  <- iemseqindex <- iemindex  <- oemvrindex <- sagaindex  <- matrix(0,nrow=length(datasizes),ncol=nsim)

for (i in (1:length(datasizes))){
  n <- datasizes[i]
  print(n)
  K <- n*20

  weight<-c(0.2, 0.8)
  mean <- 0.5
  mu<-c(mean,-mean)
  sigma<-c(1,1)*1


  weight0<-weight
  mean0 <- 1
  mu0<-c(mean0,-mean0)
  sigma0<-sigma
  seed0=23422

  ylim <- c(0.3)

  M <- 1

  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  theta<-list(p=weight,mu=mu,sigma=sigma)
  theta0<-list(p=weight0,mu=mu0,sigma=sigma0)
  # theta0<-theta

  x <- matrix(0,nrow=n,ncol=nsim)
  mls <- list()
  end <- 200
  for (j in (1:nsim))
  {

    print(j)
    seed <- j*seed0
    set.seed(seed)
    xj<-mixt.simulate(n,weight,mu,sigma)
    x[,j] <- xj
    

    df <- mixt.em(x[,j], theta0, end)
    a1 = c(rep(df[end,2],(K+1)))
    a2 = c(rep(df[end,3],(K+1)))
    b1 = c(rep(df[end,4],(K+1)))
    b2 = c(rep(df[end,5],(K+1)))
    d1 = c(rep(df[end,6],(K+1)))
    d2 = c(rep(df[end,7],(K+1)))

    ML <- cbind(1:(K+1),a1,a2,b1,b2,d1,d2)
    mls[[j]] <- ML
  }


  print('EM')
  dem <- NULL

  df.em <- vector("list", length=nsim)
  Kem <- K/n

  nbr<-1
  diem <- NULL
  df.iem <- vector("list", length=nsim)

  diemseq <- NULL
  df.iemseq <- vector("list", length=nsim)


  doemvr <- NULL
  df.oemvr <- vector("list", length=nsim)

  dsaga <- NULL
  df.saga <- vector("list", length=nsim)



  # rho.oemvr <-0.03
  # rho.saga <- 0.0016

  rho.oemvr <- 1/n**(2/3)
  rho.saga <- 1/n**(2/3)


  kiter = 1:K
  precision = 1e-3
  tmpemindex  <- tmpiemseqindex <- tmpiemindex  <- c()
  tmpoemvrindex <- tmpsagaindex  <- c()
  for (j in (1:nsim))
  {

    print(j)
    seed <- j*seed0
    set.seed(seed)
    ML <- mls[[j]]
    print("ML calculation done")

    df <- mixt.em(x[,j], theta0, Kem)
    # ML <- df
    # ML[1:(K+1),2:7]<- df[(K+1),2:7]
    df[,2:7] <- (df[,2:7] - ML[1:(Kem+1),2:7])^2
    tmpemindex[j] <- datasizes[[i]]*which(df[,c(4)] < precision)[1]
    df$rep <- j
    dem <- rbind(dem,df)
    df$rep <- NULL
    df.em[[j]] <- df
    print('em done')

    df <- mixt.iem(x[,j], theta0, K,nbr)
    df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    tmpiemindex[j] <- which(df[,c(4)] < precision)[1]
    df$rep <- j
    diem <- rbind(diem,df)
    df$rep <- NULL
    df.iem[[j]] <- df
    print('iem done')

    df <- mixt.iem.seq(x[,j], theta0, K,nbr)
    df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    tmpiemseqindex[j] <- which(df[,c(4)] < precision)[1]
    df$rep <- j
    diemseq <- rbind(diemseq,df)
    df$rep <- NULL
    df.iemseq[[j]] <- df
    print('iemseq done')

   
    df <- mixt.oemvr(x[,j], theta0, K,nbr,rho.oemvr)
    df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    tmpoemvrindex[j] <- which(df[,c(4)] < precision)[1]
    df$rep <- j
    doemvr <- rbind(doemvr,df)
    df$rep <- NULL
    df.oemvr[[j]] <- df
    print('oemvr done')

    df <- mixt.saga(x[,j], theta0, K,nbr,rho.saga)
    df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    tmpsagaindex[j] <- which(df[,c(4)] < precision)[1]
    df$rep <- j
    dsaga <- rbind(dsaga,df)
    df$rep <- NULL
    df.saga[[j]] <- df
    print('saga done')

  }

  emindex[i,] = tmpemindex
  iemindex[i,] = tmpiemindex
  iemseqindex[i,] = tmpiemseqindex
  oemvrindex[i,] = tmpoemvrindex
  sagaindex[i,] = tmpsagaindex

}




finalem <- data.frame(datasize=datasizes, emindex,algo='EM')
finalem <- melt(finalem, id.var = c("datasize","algo"))
finaliem <- data.frame(datasize=datasizes, iemindex,algo='IEM')
finaliem <- melt(finaliem, id.var = c("datasize","algo"))
finaliemseq <- data.frame(datasize=datasizes, iemseqindex,algo='IEMseq')
finaliemseq <- melt(finaliemseq, id.var = c("datasize","algo"))
finaloemvr <- data.frame(datasize=datasizes, oemvrindex,algo='OEMVR')
finaloemvr <- melt(finaloemvr, id.var = c("datasize","algo"))
finalsaga <- data.frame(datasize=datasizes, sagaindex,algo='SAGA')
finalsaga <- melt(finalsaga, id.var = c("datasize","algo"))
colnames(finalem) <-colnames(finaliemseq) <- colnames(finaliem) <- c("datasize","algo","variable","iter")
colnames(finaloemvr) <-colnames(finalsaga) <- c("datasize","algo","variable","iter")

df <- rbind(finalem,finaliem,finaliemseq,finaloemvr,finalsaga)

write.csv(df, file = "precisionsim2_moreindepruns.csv")

