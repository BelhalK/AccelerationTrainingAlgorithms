require(ggplot2)
require(gridExtra)
require(reshape2)
library(rlist)

source("utils/algos.R")
source("utils/func.R")
source("utils/plots.R")
theme_set(theme_bw())
options(digits = 2)

n <- 100000
K <- n*10
nsim=1

weight<-c(0.2, 0.8)
mean <- 0.5
mu<-c(mean,-mean)
sigma<-c(1,1)*1


weight0<-weight
mean0 <- 0.4
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

doem <- NULL
df.oem <- vector("list", length=nsim)


dsaga <- NULL
df.saga <- vector("list", length=nsim)


doemvr <- NULL
df.oemvr <- vector("list", length=nsim)

rho.oemvr <- 0.003
rho.saga <-  0.003

# rho.oemvr <- 1/n**(2/3)

# rho.saga <- 1/n**(2/3)

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
  df[,2:7] <- (df[,2:7] - ML[1:(Kem+1),2:7])^2
  # df[,2:7] <- (df[,2:7] - ML[,2:7])^2
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
iemseq$algo <- 'IEM seq'
oem$algo <- 'OEM'
oemvr$algo <- 'SVR-EM'
saga$algo <- 'FI-EM'
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
oem_ep <- oem[epochs,]
oem_ep$iteration <- 1:(K/n)
oemvr_ep <- oemvr[epochs,]
oemvr_ep$iteration <- 1:(K/n)
saga_ep <- saga[epochs,]
saga_ep$iteration <- 1:(K/n)


oemvr_ep

em_ep
# df <- NULL
# em <- iem <- oem <- oemvr <- saga <- NULL
# dem <- NULL

# df.em <- vector("list", length=nsim)
# Kem <- K/n

# nbr<-1
# diem <- NULL
# df.iem <- vector("list", length=nsim)

# doem <- NULL
# df.oem <- vector("list", length=nsim)

# doemvr <- NULL
# df.oemvr <- vector("list", length=nsim)

# dsaga <- NULL
# df.saga <- vector("list", length=nsim)

# save.image("RData/test_bign.RData")


start =10000
end = 100000




em_ep <- em
em_ep$iteration <- n*em$iteration



variance <- rbind(oemvr[start:end,c(1,4,8)],
                  iem[start:end,c(1,4,8)],
                  oem[start:end,c(1,4,8)],
                  iemseq[start:end,c(1,4,8)],
                  em_ep[2:length(epochs),c(1,4,8)],saga[start:end,c(1,4,8)])

graphConvMC2_new(variance, title="",legend=TRUE)


testoemvr <- oemvr
testiem <- iem
testoem <- oem
testiemseq <- iemseq
testem_ep <- em_ep
testsaga <- saga

testoemvr$iteration <- testoemvr$iteration/n
testiem$iteration <- testiem$iteration/n
testoem$iteration <- testoem$iteration/n
testiemseq$iteration <- testiemseq$iteration/n
testsaga$iteration <- testsaga$iteration/n
testem_ep$iteration <- testem_ep$iteration/n
testiemseq$algo <- 'IEM'

start =10000
end = 100000
variance <- rbind(testoemvr[start:end,c(1,4,8)],
                  testoem[start:end,c(1,4,8)],
                  testiemseq[start:end,c(1,4,8)],
                  testem_ep[2:11,c(1,4,8)],testsaga[start:end,c(1,4,8)])

save <- graphConvMC2_new(variance, title="",legend=TRUE)

testoemvr[start:(start+15),c(1,4,8)]
write.csv(variance, file = "notebooks/singlerun.csv")

ggsave(save,file="/Users/karimimohammedbelhal/Desktop/test.pdf", width = 200, height = 200, units = "mm")

