require(ggplot2)
require(gridExtra)
require(reshape2)
library(rlist)

source("utils/algos.R")
source("utils/func.R")
# source("utils/plots.R")
# theme_set(theme_bw())
options(digits = 22)


eml <- ieml <- iemseql <- oeml <- oemvrl <- sagal <- list()
emiter <- iemiter <-iemseqiter <- oemiter <- oemvriter <- sagaiter <- list()

datasizes <- seq(1000, 10000, 1000)

nsim=5

for (i in (1:length(datasizes))){
  n <- datasizes[i]
  print(n)
  K <- n*30

  weight<-c(0.2, 0.8)
  mean <- 0.5
  mu<-c(mean,-mean)
  sigma<-c(1,1)*1


  weight0<-weight
  mean0 <- 1
  mu0<-c(mean0,-mean0)
  sigma0<-sigma
  seed0=23422


  # ylim <- c(0.15, 0.5, 0.4)
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

    # print(j)
    seed <- 1*seed0
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


  dem <- NULL

  df.em <- vector("list", length=nsim)
  Kem <- K/n

  nbr<-1
  diem <- NULL
  df.iem <- vector("list", length=nsim)

  diemseq <- NULL
  df.iemseq <- vector("list", length=nsim)

  doem <- NULL
  df.oem <- vector("list", length=nsim)

  doemvr <- NULL
  df.oemvr <- vector("list", length=nsim)

  dsaga <- NULL
  df.saga <- vector("list", length=nsim)



  # rho.oemvr <-0.03
  # rho.saga <- 0.0016

  rho.oemvr <- 1/n**(2/3)
  rho.saga <- 1/n**(2/3)


  kiter = 1:K
  rho.oem = 1/(kiter+5)


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
    df$rep <- j
    dem <- rbind(dem,df)
    df$rep <- NULL
    df.em[[j]] <- df
    print('em done')

    df <- mixt.iem(x[,j], theta0, K,nbr)
    df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    df$rep <- j
    diem <- rbind(diem,df)
    df$rep <- NULL
    df.iem[[j]] <- df
    print('iem done')

    df <- mixt.iem.seq(x[,j], theta0, K,nbr)
    df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    df$rep <- j
    diemseq <- rbind(diemseq,df)
    df$rep <- NULL
    df.iemseq[[j]] <- df
    print('iemseq done')

    df <- mixt.oem(x[,j], theta0, K,nbr,rho.oem)
    df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    df$rep <- j
    doem <- rbind(doem,df)
    df$rep <- NULL
    df.oem[[j]] <- df
    print('oem done')

    df <- mixt.oemvr(x[,j], theta0, K,nbr,rho.oemvr)
    df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    df$rep <- j
    doemvr <- rbind(doemvr,df)
    df$rep <- NULL
    df.oemvr[[j]] <- df
    print('oemvr done')

    df <- mixt.saga(x[,j], theta0, K,nbr,rho.saga)
    df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    df$rep <- j
    dsaga <- rbind(dsaga,df)
    df$rep <- NULL
    df.saga[[j]] <- df
    print('saga done')

  }


  # dem[,2:7] <- dem[,2:7]^2
  em <- NULL
  em <- dem[dem$rep==1,]

  if (nsim>2) {
     for (j in (2:nsim))
  	{
  	  em[,2:7] <- em[,2:7]+dem[dem$rep==j,2:7]
  	}
  }
  em[,2:7] <- 1/nsim*em[,2:7]
  em[,9]<-NULL



  iem <- NULL
  iem <- diem[diem$rep==1,]

  if (nsim>2) {
  		for (j in (2:nsim))
  	{
  	  iem[,2:7] <- iem[,2:7]+diem[diem$rep==j,2:7]
  	}
  }

  iem[,2:7] <- 1/nsim*iem[,2:7]
  iem[,9]<-NULL


  iemseq <- NULL
  iemseq <- diemseq[diemseq$rep==1,]

  if (nsim>2) {
      for (j in (2:nsim))
    {
      iemseq[,2:7] <- iemseq[,2:7]+diemseq[diemseq$rep==j,2:7]
    }
  }

  iemseq[,2:7] <- 1/nsim*iemseq[,2:7]
  iemseq[,9]<-NULL


  oem <- NULL
  oem <- doem[doem$rep==1,]

  if (nsim>2) {
      for (j in (2:nsim))
    {
      oem[,2:7] <- oem[,2:7]+doem[doem$rep==j,2:7]
    }
  }

  oem[,2:7] <- 1/nsim*oem[,2:7]
  oem[,9]<-NULL




  oemvr <- NULL
  oemvr <- doemvr[doemvr$rep==1,]

  if (nsim>2) {
      for (j in (2:nsim))
    {
      oemvr[,2:7] <- oemvr[,2:7]+doemvr[doemvr$rep==j,2:7]
    }
  }

  oemvr[,2:7] <- 1/nsim*oemvr[,2:7]
  oemvr[,9]<-NULL



  saga <- NULL
  saga <- dsaga[dsaga$rep==1,]

  if (nsim>2) {
      for (j in (2:nsim))
    {
      saga[,2:7] <- saga[,2:7]+dsaga[dsaga$rep==j,2:7]
    }
  }

  saga[,2:7] <- 1/nsim*saga[,2:7]
  saga[,9]<-NULL


  iem$algo <- 'IEM'
  iemseq$algo <- 'IEMseq'
  oem$algo <- 'OEM'
  oemvr$algo <- 'OEMvr'
  saga$algo <- 'saga'
  em$algo <- 'EM'

  em$rep <- NULL
  iem$rep <- NULL
  iemseq$rep <- NULL
  oem$rep <- NULL
  oemvr$rep <- NULL
  saga$rep <- NULL


  ### PER EPOCH
  epochs = seq(1, K, by=n)
  em_ep <- em[1:(K/n),]
  em_ep$iteration <- 1:(K/n)
  iem_ep <- iem[epochs,]
  iem_ep$iteration <- 1:(K/n)
  iemseq_ep <- iemseq[epochs,]
  iemseq_ep$iteration <- 1:(K/n)
  oem_ep <- oem[epochs,]
  oem_ep$iteration <- 1:(K/n)
  oemvr_ep <- oemvr[epochs,]
  oemvr_ep$iteration <- 1:(K/n)
  saga_ep <- saga[epochs,]
  saga_ep$iteration <- 1:(K/n)


  emiter[[i]] <- em
  iemiter[[i]] <- iem
  iemseqiter[[i]] <- iemseq
  oemiter[[i]] <- oem
  oemvriter[[i]] <- oemvr
  sagaiter[[i]] <- saga

  eml[[i]] <- em_ep
  ieml[[i]] <- iem_ep
  iemseql[[i]] <- iemseq_ep
  oeml[[i]] <- oem_ep
  oemvrl[[i]] <- oemvr_ep
  sagal[[i]] <- saga_ep

}

emitersmall = emiter
iemitersmall = iemiter
iemseqitersmall = iemseqiter
oemitersmall = oemiter
oemvritersmall = oemvriter
sagaitersmall = sagaiter

# save(emitersmall,
#   iemitersmall,
#   iemseqitersmall,
#   oemitersmall,
#   oemvritersmall,
#   sagaitersmall, file = "precisionagainstn_smalln.RData")



eqiem = function(x){x*5}
eqsaga = function(x){(x*10)**(2/3)}

precision = 1e-3
emindex  <- iemseqindex <- iemindex  <- oemvrindex <- sagaindex  <- c()

for (i in (1:length(datasizes))){
  emindex[i] <- datasizes[[i]]*which(emiter[[i]][,c(4)] < precision)[1]
  iemindex[i] <- which(iemiter[[i]][,c(4)] < precision)[1]
  iemseqindex[i] <- which(iemseqiter[[i]][,c(4)] < precision)[1]
  sagaindex[i] <- which(sagaiter[[i]][,c(4)] < precision)[1]
  oemvrindex[i] <- which(oemvriter[[i]][,c(4)] < precision)[1]
}

x  <- datasizes
y1 <- data.frame(iteration=datasizes, prec = iemindex, algo = "IEM")
y2 <- data.frame(iteration=datasizes, prec = iemseqindex, algo = "IEM")
y3 <- data.frame(iteration=datasizes, prec = sagaindex, algo = "FI-EM")
y4 <- data.frame(iteration=datasizes, prec = oemvrindex, algo = "SVR-EM")
y5 <- data.frame(iteration=datasizes, prec = emindex, algo = "EM")

ytwothird <- data.frame(iteration=datasizes, prec = eqsaga(datasizes), algo = "f(n) = n^(2/3)")
ylin <- data.frame(iteration=datasizes, prec = eqiem(datasizes), algo = "f(n) = n")

sublin = eqsaga(datasizes)
lin = eqiem(datasizes)
list(sublin)[[1]]

df <- rbind(y2,y3,y4,y5,ytwothird,ylin)
colnames(df) <- c("iteration","prec","algo")


plotn <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df,aes(colour=df$algo ,linetype=df$algo))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend) + 
    scale_linetype_manual(values = c("solid","solid","solid","solid","dotted","dotted")) +
      xlab("epochs")+ scale_y_log10()  + scale_x_log10() +xlab("Problem size n") + ylab("")+
      guides(color = guide_legend(override.aes = list(size = 2))) +
      theme(axis.text.x = element_text(face="bold", color="black", 
                           size=20, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=20, angle=0),axis.title = element_text( color="black", face="bold",size=20),legend.text=element_text(size=20),
          legend.title = element_blank(),legend.position = c(0.2, 0.8))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=1))
}

plotn(df)


save <- plotn(df)
