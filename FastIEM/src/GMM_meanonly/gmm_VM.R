library(rlist)
source("utils/algos.R")
source("utils/func.R")
options(digits = 22)


n <- 100000
K <- n*100
nsim=1

weight<-c(0.2, 0.8)
mu<-c(1,-1)
sigma<-c(1,1)*1


weight0<-weight
mu0<-c(0.9,-0.9)
sigma0<-sigma
seed0=44444


# ylim <- c(0.15, 0.5, 0.4)
ylim <- c(0.3)

M <- 1
#
G<-length(mu)
col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
theta<-list(p=weight,mu=mu,sigma=sigma)
theta0<-list(p=weight0,mu=mu0,sigma=sigma0)
# theta0<-theta

a1 = c(rep(theta$p[1],(K+1)))
a2 = c(rep(theta$p[2],(K+1)))
b1 = c(rep(theta$mu[1],(K+1)))
b2 = c(rep(theta$mu[2],(K+1)))
d1 = c(rep(theta$sigma[1],(K+1)))
d2 = c(rep(theta$sigma[2],(K+1)))

ML <- cbind(1:(K+1),a1,a2,b1,b2,d1,d2)

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

rho.oemvr <- 1/n**(2/3)
rho.saga <- 1/n**(2/3)
kiter = 1:K
rho.oem = 3/(kiter+10)


x <- matrix(0,nrow=n,ncol=nsim)


for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  xj<-mixt.simulate(n,weight,mu,sigma)
  x[,j] <- xj

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
oem_ep <- oem[epochs,]
oem_ep$iteration <- 1:(K/n)
oemvr_ep <- oemvr[epochs,]
oemvr_ep$iteration <- 1:(K/n)
saga_ep <- saga[epochs,]
saga_ep$iteration <- 1:(K/n)



# start =0
# end = 10
# variance <- rbind(oemvr_ep[start:end,c(1,5,8)],iem_ep[start:end,c(1,5,8)],
#                   oem_ep[start:end,c(1,5,8)],em_ep[start:end,c(1,5,8)],saga_ep[start:end,c(1,5,8)])



df <- NULL
em <- iem <- oem <- oemvr <- saga <- NULL
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

save.image("RData/gmm_big2.RData")