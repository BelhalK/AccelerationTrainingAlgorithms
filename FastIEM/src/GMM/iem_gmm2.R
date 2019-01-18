source("algos.R")
source("func.R")
source("plots.R")
theme_set(theme_bw())
options(digits = 22)
n <- 100
weight<-c(0.4, 0.6) 
mu<-c(-2,2)
sigma<-c(1,1)*1


weight0<-c(.5,.5)
mu0<-c(-5,5)
sigma0<-c(0.8,0.8)


K <- 1000

seed0=44444


# ylim <- c(0.15, 0.5, 0.4)
ylim <- c(0.3)

M <- 1
nsim <- 3
#
G<-length(mu)
col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
theta<-list(p=weight,mu=mu,sigma=sigma)
theta0<-list(p=weight0,mu=mu0,sigma=sigma0)
# theta0<-theta


##  Simulation
x <- matrix(0,nrow=n,ncol=nsim)
for (j in (1:nsim))
{
  seed <- j*seed0
  set.seed(seed)
  xj<-mixt.simulate(n,weight,mu,sigma)
  x[,j] <- xj
}


## EM
print('EM')
dem <- NULL
df.em <- vector("list", length=nsim)
Kem <- K/n
for (j in (1:nsim))
{ print(j)
  df <- mixt.em(x[,j], theta0, K)
  df <- mixt.ident(df)
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}
graphConvMC_new(dem, title="EM")



## EM
print('EM')
diem <- NULL
df.iem <- vector("list", length=nsim)
for (j in (1:nsim))
{ print(j)
  df <- mixt.iem(x[,j], theta0, K,1)
  df <- mixt.ident(df)
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}
graphConvMC_new(diem, title="EM")


## EM
print('EM')
doemvr <- NULL
df.oemvr <- vector("list", length=nsim)
nbr <- 1
for (j in (1:nsim))
{ print(j)
  df <- mixt.oemvr(x[,j], theta0, K,nbr,0.003)
  df <- mixt.ident(df)
  df$rep <- j
  doemvr <- rbind(doemvr,df)
  df$rep <- NULL
  df.oemvr[[j]] <- df
}
graphConvMC_new(doemvr, title="EM")

################################################

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
rho <- 0.0003
for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  x <- matrix(0,nrow=n,ncol=nsim)
  xj<-mixt.simulate(n,weight,mu,sigma)
  x[,j] <- xj

  df <- mixt.em(x[,j], theta0, Kem)
  # ML <- df
  # ML[1:(K+1),2]<- df[(K+1),2]
  df[1:Kem,2:7] <- (df[1:Kem,2:7] - ML[1:Kem,2:7])^2
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df

  df <- mixt.iem(x[,j], theta0, K,nbr)
  df[,2:7] <- (df[,2:7] - ML[,2:7])^2
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df

  df <- mixt.oem(x[,j], theta0, K,nbr)
  df[,2:7] <- (df[,2:7] - ML[,2:7])^2
  df$rep <- j
  doem <- rbind(doem,df)
  df$rep <- NULL
  df.oem[[j]] <- df

  df <- mixt.oemvr(x[,j], theta0, K,nbr,rho)
  df[,2:7] <- (df[,2:7] - ML[,2:7])^2
  df$rep <- j
  doemvr <- rbind(doemvr,df)
  df$rep <- NULL
  df.oemvr[[j]] <- df
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

iem$algo <- 'IEM'
oem$algo <- 'OEM'
oemvr$algo <- 'OEMvr'
em$algo <- 'EM'

em$rep <- NULL
iem$rep <- NULL
oem$rep <- NULL
oemvr$rep <- NULL




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
# m <- graphConvMC(em_scaled[0:K,c(1,2,8)],iem[0:K,c(1,2,8)],oem[0:K,c(1,2,8)])
# m <- graphConvMC(em_scaled[0:K,c(1,2,9)],iem[0:K,c(1,2,9)],oem[0:K,c(1,2,9)],oemvr[0:K,c(1,2,9)])

variance <- NULL
variance <- rbind(em_ep[1:(K/n),],iem_ep[1:(K/n),],oemvr_ep[1:(K/n),], oem_ep[1:(K/n),])
# variance <- rbind(iem_ep[2:(K/n),],oemvr_ep[2:(K/n),])
variance <- rbind(em_ep[2:(K/n),],iem_ep[2:(K/n),],oemvr_ep[2:(K/n),], oem_ep[2:(K/n),])
graphConvMC2_new(variance, title="IEMs",legend=TRUE)
