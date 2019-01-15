source("algos.R")
source("func.R")
source("plots.R")
theme_set(theme_bw())
options(digits = 22)
# save.image("gmm_mu.RData")
load("gmm_mu.RData")
n <- 1000
weight<-c(0.2, 0.8) 
mu<-c(1,-1)
sigma<-c(1,1)*1


weight0<-weight
mu0<-c(2,-0.5)
sigma0<-sigma


K <- 2000

seed0=44444


# ylim <- c(0.15, 0.5, 0.4)
ylim <- c(0.3)

M <- 1
nsim <- 5
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


# ## EM
# print('EM')
# doemvr <- NULL
# df.oemvr <- vector("list", length=nsim)

# for (j in (1:nsim))
# { print(j)
#   df <- mixt.oemvr(x[,j], theta0, K,nbr)
#   df <- mixt.ident(df)
#   df$rep <- j
#   doemvr <- rbind(doemvr,df)
#   df$rep <- NULL
#   df.oemvr[[j]] <- df
# }
# graphConvMC_new(doemvr, title="EM")

################################################
ML <- df
for (i in (1:(K+1))){
  ML[i,2:7]<- c(theta$p,theta$mu,theta$sigma)
}

print('EM')
dem <- NULL

df.em <- vector("list", length=nsim)

nbr<-1
diem <- NULL
df.iem <- vector("list", length=nsim)

doem <- NULL
df.oem <- vector("list", length=nsim)

doemvr <- NULL
df.oemvr <- vector("list", length=nsim)

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  x <- matrix(0,nrow=n,ncol=nsim)
  xj<-mixt.simulate(n,weight,mu,sigma)
  x[,j] <- xj

  df <- mixt.em(x[,j], theta0, K)
  # ML <- df
  # ML[1:(K+1),2:7]<- df[(K+1),2:7]
  df[,2:7] <- (df[,2:7] - ML[,2:7])^2
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

  df <- mixt.oemvr(x[,j], theta0, K,nbr)
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


em_scaled <- em
em_scaled$iteration = seq(0, n*K, by=n)
em_scaled <- em_scaled[rep(seq_len(nrow(em_scaled)), each=n),]


iem$algo <- 'IEM'
oem$algo <- 'OEM'
oemvr$algo <- 'OEMvr'
em_scaled$algo <- 'EM'

em_scaled$rep <- NULL
iem$rep <- NULL
oem$rep <- NULL
oemvr$rep <- NULL


variance <- NULL
variance <- rbind(oemvr[1001:2001,c(1,5,8)],iem[1001:2001,c(1,5,8)],oem[1001:2001,c(1,5,8)],em_scaled[1001:2001,c(1,5,8)])
variance <- rbind(oemvr[2:(K+1),c(1,5,8)],iem[2:(K+1),c(1,5,8)],oem[2:(K+1),c(1,5,8)],em_scaled[2:(K+1),c(1,5,8)])
graphConvMC2_new(variance, title="IEMs",legend=TRUE)


