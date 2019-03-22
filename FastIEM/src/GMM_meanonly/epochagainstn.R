require(ggplot2)
require(gridExtra)
require(reshape2)
library(rlist)

source("utils/algos.R")
source("utils/func.R")
source("utils/plots.R")
theme_set(theme_bw())
options(digits = 22)

# save.image("RData/precisionagainstn.RData")
# load("RData/precisionagainstn_VM.RData")
load("RData_VM/precisionagainstn_VM.RData")

length(datasizes)
eml <- ieml <- oeml <- oemvrl <- sagal <- list()
datasizes <- c(1000, 5000, 10000, 20000, 30000, 40000,50000)

nsim=1

for (i in (1:length(datasizes))){
  n <- datasizes[i]
  K <- n*10

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

    print(j)
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


  print('EM')
  dem <- NULL

  df.em <- vector("list", length=nsim)
  Kem <- K/n

  nbr<-1
  diem <- NULL
  df.iem <- vector("list", length=nsim)

  doem <- NULL
  df.oem <- vector("list", length=nsim)

  doemvr <- NULL
  df.oemvr <- vector("list", length=nsim)

  dsaga <- NULL
  df.saga <- vector("list", length=nsim)



  rho.oemvr <-0.03
  # rho.saga <- 0.0016

  # rho.oemvr <- 1/n**(2/3)
  rho.saga <- 1/n**(2/3)
  # rho.saga <- 0.001

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

    # df <- mixt.oem(x[,j], theta0, K,nbr,rho.oem)
    # df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    # df$rep <- j
    # doem <- rbind(doem,df)
    # df$rep <- NULL
    # df.oem[[j]] <- df
    # print('oem done')

    # df <- mixt.oemvr(x[,j], theta0, K,nbr,rho.oemvr)
    # df[,2:7] <- (df[,2:7] - ML[,2:7])^2
    # df$rep <- j
    # doemvr <- rbind(doemvr,df)
    # df$rep <- NULL
    # df.oemvr[[j]] <- df
    # print('oemvr done')

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
  oem$algo <- 'OEM'
  oemvr$algo <- 'OEMvr'
  saga$algo <- 'saga'
  em$algo <- 'EM'

  em$rep <- NULL
  iem$rep <- NULL
  oem$rep <- NULL
  oemvr$rep <- NULL
  saga$rep <- NULL


  ### PER EPOCH
  epochs = seq(1, K, by=n)
  em_ep <- em[1:(K/n),]
  em_ep$iteration <- 1:(K/n)
  iem_ep <- iem[epochs,]
  iem_ep$iteration <- 1:(K/n)
  # oem_ep <- oem[epochs,]
  # oem_ep$iteration <- 1:(K/n)
  # oemvr_ep <- oemvr[epochs,]
  # oemvr_ep$iteration <- 1:(K/n)
  saga_ep <- saga[epochs,]
  saga_ep$iteration <- 1:(K/n)

  eml[[i]] <- em_ep
  ieml[[i]] <- iem_ep
  # oeml[[i]] <- oem_ep
  # oemvrl[[i]] <- oemvr_ep
  sagal[[i]] <- saga_ep

}


# precision = 0.00001
# emindex <- iemindex <- oemindex <- oemvrindex <- sagaindex <- c()

# for (i in (1:length(datasizes))){
#   emindex[i] <- which(eml[[i]][,c(4)] < precision)[1]
#   iemindex[i] <- which(ieml[[i]][,c(4)] < precision)[1]
#   # oemindex[i] <- which(oeml[[i]][,c(4)] < precision)[1]
#   # oemvrindex[i] <- which(oemvrl[[i]][,c(4)] < precision)[1]
#   sagaindex[i] <- which(sagal[[i]][,c(4)] < precision)[1]
# }


# emindex
# iemindex
# oemindex
# oemvrindex
# sagaindex


precision = 1e-7
for (precision in c(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11)){
  emindex <- iemindex  <- sagaindex <- c()

  for (i in (1:length(datasizes))){
    emindex[i] <- datasizes[[i]]*which(emiter[[i]][,c(4)] < precision)[1]
    iemindex[i] <- which(iemiter[[i]][,c(4)] < precision)[1]
    sagaindex[i] <- which(sagaiter[[i]][,c(4)] < precision)[1]
  }
  print(precision)
  print("EM")
  print(emindex)
  print("IEM")
  print(iemindex)
  print("SAGA")
  print(sagaindex)

  x  <- datasizes
  y1 <- emindex
  y2 <- iemindex
  y3 <- sagaindex
  df <- data.frame(x,y1,y2, y3)

  print(ggplot(df, aes(x),show.legend = TRUE) +                    
  geom_line(aes(y=y1), colour="red") +  
  geom_line(aes(y=y2), colour="green") +
  geom_line(aes(y=y3), colour="purple") +
  xlab("Dataset size") + ylab("Epoch")  +
  ggtitle(precision))

}



eqiem = function(x){x}
eqsaga = function(x){x**(2/3)}

x  <- datasizes
y1 <- eqiem(datasizes)
y2 <- eqsaga(datasizes)
df <- data.frame(x,y1,y2)

ggplot(df, aes(x),show.legend = TRUE) +                    
  geom_line(aes(y=y1), colour="red") +  
  geom_line(aes(y=y2), colour="green") +
  xlab("Dataset size") + ylab("Epoch")  +
  ggtitle(precision)


plot(eqsaga(1000:100000), type='l')



for (precision in c(1e-2,1e-3,1e-4,1e-5)){
  emindex <- iemindex  <- sagaindex <- c()

  for (i in (1:length(datasizes))){
    iemindex[i] <- which(iemiter[[i]][,c(4)] < precision)[1]
    iemseqindex[i] <- which(iemseqiter[[i]][,c(4)] < precision)[1]
    sagaindex[i] <- which(sagaiter[[i]][,c(4)] < precision)[1]
    oemvrindex[i] <- which(oemvriter[[i]][,c(4)] < precision)[1]
  }
  print("IEM")
  print(iemindex)
  print("SAGA")
  print(sagaindex)

  x  <- datasizes
  y1 <- iemindex
  y2 <- iemseqindex
  y3 <- sagaindex
  y4 <- oemvrindex
  ytwothird <- eqsaga(datasizes)
  ylin <- eqiem(datasizes)
  df <- data.frame(x,y2, y3,ytest)

  print(ggplot(df, aes(x),show.legend = TRUE) +                    
    geom_line(aes(y=y1), colour="red") +
  geom_line(aes(y=y2), colour="blue") +
  geom_line(aes(y=y3), colour="green") +
  geom_line(aes(y=y4), colour="purple") +
  geom_line(aes(y=ytwothird), colour="black") +
  geom_line(aes(y=ylin), colour="brown") +
  scale_y_log10()  + scale_x_log10() +
  xlab("Dataset size") + ylab("Epoch")  +
  ggtitle(precision))

}




for (precision in c(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11)){
  emindex <- iemindex  <- sagaindex <- c()

  for (i in (1:length(datasizes))){
    sagaindex[i] <- which(sagaiter[[i]][,c(4)] < precision)[1]
  }

  print("SAGA")
  print(sagaindex)

  x  <- datasizes
  y3 <- sagaindex
  df <- data.frame(x,y2, y3)

  print(ggplot(df, aes(x),show.legend = TRUE) +                    
  geom_line(aes(y=y3), colour="purple") +
  xlab("Dataset size") + ylab("Epoch")  +
  ggtitle(precision))

}


curvesaga <- sagal
curveem <- eml
curveiem <- ieml

for (i in (1:length(datasizes))){
  curvesaga[[i]]$size <- datasizes[i]
  curveem[[i]]$size <- datasizes[i]
  curveiem[[i]]$size <- datasizes[i]
}



epochs
start =1
end = 10

plot.saga <- NULL
for (i in (1:length(datasizes))){
  plot.saga <- rbind(plot.saga, curvesaga[[i]][start:end,c(1,4,9)])

}

# plot.saga <- rbind(curvesaga[[1]][start:end,c(1,4,9)],
#                   curvesaga[[2]][start:end,c(1,4,9)],
#                   curvesaga[[3]][start:end,c(1,4,9)],
#                   curvesaga[[4]][start:end,c(1,4,9)],
#                   curvesaga[[5]][start:end,c(1,4,9)],
#                   curvesaga[[6]][start:end,c(1,4,9)],
#                   curvesaga[[7]][start:end,c(1,4,9)],
#                   curvesaga[[8]][start:end,c(1,4,9)])
plotagainstn(plot.saga, title="SAGA-EM for different dataset size",legend=TRUE)


plot.iem <- rbind(curveiem[[1]][start:end,c(1,4,9)],
                  curveiem[[2]][start:end,c(1,4,9)],
                  curveiem[[3]][start:end,c(1,4,9)],
                  curveiem[[4]][start:end,c(1,4,9)],
                  curveiem[[5]][start:end,c(1,4,9)],
                  curveiem[[6]][start:end,c(1,4,9)],
                  curveiem[[7]][start:end,c(1,4,9)],
                  curveiem[[8]][start:end,c(1,4,9)])
plotagainstn(plot.iem, title="IEMs GMM 1e5",legend=TRUE)



plot.em <- rbind(curveem[[1]][start:end,c(1,4,9)],
                  curveem[[2]][start:end,c(1,4,9)],
                  curveem[[3]][start:end,c(1,4,9)],
                  curveem[[4]][start:end,c(1,4,9)],
                  curveem[[5]][start:end,c(1,4,9)],
                  curveem[[6]][start:end,c(1,4,9)],
                  curveem[[7]][start:end,c(1,4,9)],
                  curveem[[8]][start:end,c(1,4,9)])
plotagainstn(plot.em, title="IEMs GMM 1e5",legend=TRUE)



epochs
start =1
end = 20

i = 4
variance <- rbind(ieml[[i]][start:end,c(1,4,8)],
                  eml[[i]][start:end,c(1,4,8)],
                   sagal[[i]][start:end,c(1,4,8)])

graphConvMC2_new(variance, title="IEMs GMM 1e5",legend=TRUE)


# epochs
# start =7
# end = 10
# variance <- rbind(oemvr_ep[start:end,c(1,4,8)],
#                   iem_ep[start:end,c(1,4,8)],
#                   em_ep[start:end,c(1,4,8)],
#                    saga_ep[start:end,c(1,4,8)])

# graphConvMC2_new(variance, title="IEMs",legend=TRUE)


# epochs
# start =2
# end = 10

# variance <- rbind(oemvr_ep[start:end,c(1,4,8)],
#                   iem_ep[start:end,c(1,4,8)],
#                   oem_ep[start:end,c(1,4,8)],
#                   em_ep[start:end,c(1,4,8)])


# graphConvMC2_new(variance, title="IEMs",legend=TRUE)

# graphConvMC2_new(iem_ep[start:end,c(1,4,8)], title="IEMs",legend=TRUE)
# graphConvMC2_new(oem_ep[start:end,c(1,4,8)], title="IEMs",legend=TRUE)
# graphConvMC2_new(oemvr_ep[start:end,c(1,4,8)], title="IEMs",legend=TRUE)
# graphConvMC2_new(saga_ep[start:end,c(1,4,8)], title="IEMs",legend=TRUE)
