source("alogs.R")
source("func.R")
source("plots.R")
theme_set(theme_bw())
library(MASS)

ni <- 6
n <- 1000
mu<-c(4,9)
mu0<-c(1,4)

da <- length(mu)
db <- 3

A <- matrix(sample.int(10, size = ni*ni, replace = TRUE), nrow = ni, ncol = da)
B <- matrix(sample.int(10, size = ni*ni, replace = TRUE), nrow = ni, ncol = db)
id <- rep(1:n,each=ni)
sigma <- 0.1*diag(ni)
omega <- 0.1*diag(db)


K<-100
seed0=44444
nsim <- 3

col.names <- c("iteration", paste0("mu",1:da))
theta<-list(mu=mu)
theta0<-list(mu=mu0)
##  Simulation
x <- matrix(0,nrow=n*ni,ncol=nsim)
z <- matrix(0,nrow=n,ncol=nsim)
for (j in (1:nsim))
{
  seed <- j*seed0
  set.seed(seed)
  # xj<-mixt.simulate(n,ni,A,B,mu,sigma,omega)
  xj <- mvrnorm(n, A%*%mu,B%*%omega%*%t(B)+sigma)
}


## EM
print('EM')
dem <- NULL
nbr<-n
df.em <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.em(x[,j], theta0, K,A,B,sigma,omega,id)
  # ML <- df
  # ML[1:(K+1),2:3]<- df[(K+1),2:3]
  # df[,2:3] <- df[,2:3] - ML[,2:3]
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}

# dem[,2:3] <- dem[,2:3]^2
graphConvMC_new(dem, title="EM")

em <- NULL
em <- dem[dem$rep==1,]

if (nsim>2) {
   for (j in (2:nsim))
	{
	  em[,2:3] <- em[,2:3]+dem[dem$rep==j,2:3]
	}
}
em[,2:3] <- 1/nsim*em[,2:3]
em[,4]<-NULL


## IEM
print('IEM')
nbr<-1
diem <- NULL
df.iem <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.iem(x[,j], theta0, K,A,B,sigma,omega,id,nbr)
  ML <- df
  ML[1:(K+1),2:3]<- df[(K+1),2:3]
  df[,2:3] <- df[,2:3] - ML[,2:3]
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}
diem[,2:3] <- diem[,2:3]^2
graphConvMC_new(diem, title="IEM 1R")

iem1 <- NULL
iem1 <- diem[diem$rep==1,]

if (nsim>2) {
		for (j in (2:nsim))
	{
	  iem1[,2:3] <- iem1[,2:3]+diem[diem$rep==j,2:3]
	}
}

iem1[,2:3] <- 1/nsim*iem1[,2:3]
iem1[,4]<-NULL

## IEM
print('IEM')
nbr<-n/2
diem50 <- NULL
df.iem50 <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.iem(x[,j], theta0, K,A,B,sigma,omega,id,nbr)
  ML <- df
  ML[1:(K+1),2:3]<- df[(K+1),2:3]
  df[,2:3] <- df[,2:3] - ML[,2:3]
  df$rep <- j
  diem50 <- rbind(diem50,df)
  df$rep <- NULL
  df.iem50[[j]] <- df
}
diem50[,2:3] <- diem50[,2:3]^2
# graphConvMC_new(diem50, title="IEM 1R")

iem2 <- NULL
iem2 <- diem50[diem50$rep==1,]

if (nsim>2) {
    for (j in (2:nsim))
  {
    iem2[,2:3] <- iem2[,2:3]+diem50[diem50$rep==j,2:3]
  }
}

iem2[,2:3] <- 1/nsim*iem2[,2:3]
iem2[,4]<-NULL


em_scaled <- em
em_scaled$iteration = seq(0, n*K, by=n)
em_scaled <- em_scaled[rep(seq_len(nrow(em_scaled)), each=n),]

iem2_scaled <- iem2
iem2_scaled$iteration = seq(0, n/2*K, by=n/2)
iem2_scaled <- iem2_scaled[rep(seq_len(nrow(iem2_scaled)), each=n/2),]


iem1$algo <- 'IEM 1%'
iem2_scaled$algo <- 'IEM 50%'
em_scaled$algo <- 'EM'
variance <- NULL
variance <- rbind(em_scaled[0:K,],iem2_scaled[0:K,],iem1[0:K,])
colnames(variance) <- c("iteration","beta1","beta2","algo")
# graphConvMC2_new(variance, title="IEMs alpha=0.33",legend=TRUE)
plot_new3(variance,legend=FALSE)




############ FEW RUNS ################
K=1000
nsim=2
mus0<-list(c(1,5),c(3,7),c(2,6))

## EM
print('EM')
dem <- NULL
nbr<-n
df.em <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  theta0<-list(mu=mus0[[j]])
  df <- mixt.iem(x[,j], theta0, K,A,B,sigma,omega,id,nbr)
  df$rep <- j
  df_scaled <- df
  df_scaled$iteration = seq(0, n*K, by=n)
  df_scaled <- df_scaled[rep(seq_len(nrow(df_scaled)), each=n),]
  dem <- rbind(dem,df_scaled[1:(K+1),])
  df$rep <- NULL
  df.em[[j]] <- df
}


## IEM
print('IEM')
nbr<-n/2
diem50 <- NULL
df.iem50 <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  theta0<-list(mu=mus0[[j]])
  df <- mixt.iem(x[,j], theta0, K,A,B,sigma,omega,id,nbr)
  df$rep <- j
  df_scaled <- df
  df_scaled$iteration = seq(0, n/2*K, by=n/2)
  df_scaled <- df_scaled[rep(seq_len(nrow(df_scaled)), each=n/2),]
  diem50 <- rbind(diem50,df_scaled[1:(K+1),])
  df$rep <- NULL
  df.iem50[[j]] <- df
}


## IEM
print('IEM')
nbr<-1
diem <- NULL
df.iem <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  theta0<-list(mu=mus0[[j]])
  df <- mixt.iem(x[,j], theta0, K,A,B,sigma,omega,id,nbr)
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}

colnames(dem) <- colnames(diem50)<- colnames(diem50)<- c("iteration","beta1","beta2","rep")
graphConvMCdf3_new(dem,diem50,diem)

# graphConvMCdf2_new(diem,diem50, title="IEM 1R")



# ## IEM
# print('IEM')
# nbr<-15
# KR <- KNR*n/nbr
# diem <- NULL
# df.iem <- vector("list", length=nsim)
# for (j in (1:nsim))
# {
#   print(j)
#   df <- mixt.iem(x[,j], theta0, K,A,B,sigma,omega,id,nbr)
#   df$rep <- j
#   diem <- rbind(diem,df)
#   df$rep <- NULL
#   df.iem[[j]] <- df
# }
# # graphConvMC_new(diem, title="IEM 1R")
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
# iem3[,4]<-NULL



# iem1$algo <- '1/N'
# iem2$algo <- '25%'
# iem3$algo <- '50%'
# em$algo <- 'EM'
# variance <- NULL
# end <- 30000
# variance <- rbind(em[1:end,],iem1[1:end,],iem2[1:end,],iem3[1:end,])
# variance <- rbind(em[1:end,],iem1[1:end,])
# variance <- rbind(em,iem1)
# graphConvMC2_new(variance, title="IEMs alpha=0.33",legend=TRUE)
# plot_new(variance,legend=TRUE)


K=15000
nsim=2
mus0<-list(c(1,5),c(3,7),c(2,6))


print('EM')
dem <- NULL
nbrem<-n/2
df.em <- vector("list", length=nsim)

nbriem1<-n/10
diem <- NULL
df.iem <- vector("list", length=nsim)

nbriem2<-n/5
diem50 <- NULL
df.iem50 <- vector("list", length=nsim)

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  theta0<-list(mu=mus0[[j]])
  x <- mvrnorm(n, A%*%mu,B%*%omega%*%t(B)+sigma)
  # print(head(x))
  df <- mixt.iem(x, theta0, K,A,B,sigma,omega,id,nbrem)
  df$rep <- j
  df_scaled <- df
  df_scaled$iteration = seq(0, 5*K, by=5)
  df_scaled <- df_scaled[rep(seq_len(nrow(df_scaled)), each=5),]
  dem <- rbind(dem,df_scaled[1:(K+1),])
  df$rep <- NULL
  df.em[[j]] <- df

  df <- mixt.iem(x, theta0, K,A,B,sigma,omega,id,nbriem2)
  df$rep <- j
  df_scaled <- df
  df_scaled$iteration = seq(0, 2*K, by=2)
  df_scaled <- df_scaled[rep(seq_len(nrow(df_scaled)), each=2),]
  diem50 <- rbind(diem50,df_scaled[1:(K+1),])
  df$rep <- NULL
  df.iem50[[j]] <- df

  df <- mixt.iem(x, theta0, K,A,B,sigma,omega,id,nbriem1)
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}

colnames(dem) <- colnames(diem50)<- colnames(diem50)<- c("iteration","beta1","beta2","rep")
graphConvMCdf3_new(dem,diem50,diem)
