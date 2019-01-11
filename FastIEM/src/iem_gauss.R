source("algos_scalar.R")
source("func_scalar.R")
theme_set(theme_bw())

# n <- 100
# mu<-c(1,0)
# mu0<-c(0.1,0)
# sigma<-c(0.5,0.1)*1

n <- 100
mu<-c(10,0)
mu0<-c(5,0)
# sigma<-c(10,5)*1
sigma<-c(3,10)*1

alph <- sigma[2]/(sigma[1]+sigma[2])
gamm <- 1/(1/sigma[1]+1/sigma[2])

KNR <- 100
K1<-30
# Several Chains for the same iteration
M <- 1

alpha1 <- 0.7
alpha2 <- 0.4
seed0=44444
ylim <- c(0.3)

nsim <- 2
G<-1
col.names <- c("iteration", paste0("mu",1:G))
theta<-list(mu=mu[1])
theta0<-list(mu=mu0[1])
##  Simulation
x <- matrix(0,nrow=n,ncol=nsim)
for (j in (1:nsim))
{
  seed <- j*seed0
  set.seed(seed)
  xj<-mixt.simulate(n,mu,sigma)
  x[,j] <- xj
}

## EM
print('EM')
dem <- NULL
nbr<-n
KR <- KNR*n/nbr
K<-KR
df.em <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.em(x[,j], theta0, KR, alph)
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}
graphConvMC_new(dem, title="EM")



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
  df <- mixt.iem_seq(x[,j], theta0, KR, alph,nbr)
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}
graphConvMC_new(diem, title="IEM 1R")
# diem[,2] <- diem[,2]^2
iem1 <- NULL
iem1 <- diem[diem$rep==1,]

if (nsim>2) {
		for (j in (2:nsim))
	{
	  iem1[,2] <- iem1[,2]+diem[diem$rep==j,2]
	}
}

iem1$iteration <- 0:KR
iem1[,2] <- 1/nsim*iem1[,2]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,iem1[i,2])
}


for (l in (0:(KR-1)))
{
  iem1[(l*nbr+2):((l+1)*nbr+1),2] <- Lr[l+1,]
}
iem1$iteration <- 1:(KR*nbr+1)
iem1[,3]<-NULL



## IEM
print('IEM')
nbr<-50
KR <- KNR*n/nbr
diem <- NULL
df.iem <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.iem_seq(x[,j], theta0, KR, alph,nbr)
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}
graphConvMC_new(diem, title="IEM 1R")
# diem[,2] <- diem[,2]^2
iem5 <- NULL
iem5 <- diem[diem$rep==1,]

if (nsim>2) {
    for (j in (2:nsim))
  {
    iem5[,2] <- iem5[,2]+diem[diem$rep==j,2]
  }
}

iem5$iteration <- 0:KR
iem5[,2] <- 1/nsim*iem5[,2]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,iem5[i,2])
}


for (l in (0:(KR-1)))
{
  iem5[(l*nbr+2):((l+1)*nbr+1),2] <- Lr[l+1,]
}
iem5$iteration <- 1:(KR*nbr+1)
iem5[,3]<-NULL







# ## IEM
# print('IEM')
# nbr<-1
# KR <- KNR*n/nbr
# diem <- NULL
# df.iem <- vector("list", length=nsim)
# for (j in (1:nsim))
# {
#   print(j)
#   df <- mixt.iem_perpass(x[,j], theta0, KR, alph,nbr)
#   df$rep <- j
#   diem <- rbind(diem,df)
#   df$rep <- NULL
#   df.iem[[j]] <- df
# }
# graphConvMC_new(diem, title="IEM 1R")
# # diem[,2] <- diem[,2]^2
# iem3 <- NULL
# iem3 <- diem[diem$rep==1,]

# if (nsim>2) {
#     for (j in (2:nsim))
#   {
#     iem3[,2] <- iem3[,2]+diem[diem$rep==j,2]
#   }
# }

# iem3$iteration <- 0:KR
# iem3[,2] <- 1/nsim*iem3[,2]
# Lr <- NULL
# for (i in (2:(KR+1)))
# {
#   Lr <- rbind(Lr,iem3[i,2])
# }


# for (l in (0:(KR-1)))
# {
#   iem3[(l*nbr+2):((l+1)*nbr+1),2] <- Lr[l+1,]
# }
# iem3$iteration <- 1:(KR*nbr+1)
# iem3[,3]<-NULL






iem1$algo <- '1/N'
iem2$algo <- 'rand per iter'
iem3$algo <- 'rand per pass'
iem5$algo <- '50%'
em$algo <- 'EM'
variance <- NULL
variance <- rbind(em,iem1, iem5)
variance <- rbind(em,iem1,iem2, iem3)
graphConvMC2_new(variance, title="IEMs alpha=0.33",legend=TRUE)


