source("algos.R")
source("func.R")
theme_set(theme_bw())

n <- 100
weight<-c(0.2, 0.8) 
mu<-c(10,-10)
sigma<-c(1,1)*1


weight0<-c(.5,.5)
mu0<-c(1,2)
sigma0<-c(1,1)
nb_r <- 10
KNR <- 50
K1 <-10
K <- 300

alpha1 <- 0.7
alpha2 <- 0.4
seed0=44444


# ylim <- c(0.15, 0.5, 0.4)
ylim <- c(0.1, 0.3, 0.3)

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
for (j in (1:nsim))
{ print(j)
  df <- mixt.em(x[,j], theta0, K)
  df <- mixt.ident(df)
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}
graphConvMC(dem, title="EM")



# dem[,2] <- dem[,2]^2
em <- NULL
em <- dem[dem$rep==1,]


if (nsim>2) {
   for (j in (2:nsim))
  {
    em[,2] <- em[,2]+dem[dem$rep==j,2]
  }
}

em$iteration <- 0:KNR
em[,2] <- 1/nsim*em[,2]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,em[i,2])
}


for (l in (0:(KR-1)))
{
  em[(l*nbr+2):((l+1)*nbr+1),2] <- Lr[l+1,]
}
em$iteration <- 1:(KR*nbr+1)
em[,3]<-NULL


## IEM
print('IEM')
nbr<-1
KR <- KNR*n/nbr
diem <- NULL
df.iem <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.iem(x[,j], theta0, K,nbr)
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}
graphConvMC(diem, title="IEM 1R")
# diem[,2] <- diem[,2]^2
iem <- NULL
iem <- diem[diem$rep==1,]

if (nsim>2) {
    for (j in (2:nsim))
  {
    iem[,2] <- iem[,2]+diem[diem$rep==j,2]
  }
}

iem$iteration <- 0:KR
iem[,2] <- 1/nsim*iem[,2]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,iem[i,2])
}


for (l in (0:(KR-1)))
{
  iem[(l*nbr+2):((l+1)*nbr+1),2] <- Lr[l+1,]
}
iem$iteration <- 1:(KR*nbr+1)
iem[,3]<-NULL

# oEM
print('oEM')
nbr<-1
KR <- KNR*n/nbr
diem <- NULL
df.iem <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.oem(x[,j], theta0, KR, alph,nbr)
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}
graphConvMC_new(diem, title="IEM 1R")

oem <- NULL
oem <- diem[diem$rep==1,]

if (nsim>2) {
    for (j in (2:nsim))
  {
    oem[,2] <- oem[,2]+diem[diem$rep==j,2]
  }
}

oem$iteration <- 0:KR
oem[,2] <- 1/nsim*oem[,2]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,oem[i,2])
}

for (l in (0:(KR-1)))
{
  oem[(l*nbr+2):((l+1)*nbr+1),2] <- Lr[l+1,]
}
oem$iteration <- 1:(KR*nbr+1)
oem[,3]<-NULL


# oemvr
print('oemvr')
nbr<-1
KR <- KNR*n/nbr
diem <- NULL
df.iem <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.oemvr(x[,j], theta0, KR, alph,nbr)
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}
graphConvMC_new(diem, title="IEM 1R")

oemvr <- NULL
oemvr <- diem[diem$rep==1,]

if (nsim>2) {
    for (j in (2:nsim))
  {
    oemvr[,2] <- oemvr[,2]+diem[diem$rep==j,2]
  }
}

oemvr$iteration <- 0:KR
oemvr[,2] <- 1/nsim*oemvr[,2]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,oemvr[i,2])
}

for (l in (0:(KR-1)))
{
  oemvr[(l*nbr+2):((l+1)*nbr+1),2] <- Lr[l+1,]
}
oemvr$iteration <- 1:(KR*nbr+1)
oemvr[,3]<-NULL



iem$algo <- 'iem'
oem$algo <- 'oem'
oemvr$algo <- 'oem-vr'
em$algo <- 'EM'

variance <- NULL
variance <- rbind(em,iem,oem,oemvr)
graphConvMC2_new(variance, title="IEMs alpha=0.33",legend=TRUE)