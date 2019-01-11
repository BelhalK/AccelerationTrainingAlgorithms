source("alogs.R")
source("func.R")
theme_set(theme_bw())


ni <- 6
n <- 30
mu<-c(2,8)
mu0<-c(1,5)

da <- length(mu)
db <- 3

A <- matrix(sample.int(10, size = ni*ni, replace = TRUE), nrow = ni, ncol = da)
B <- matrix(sample.int(10, size = ni*ni, replace = TRUE), nrow = ni, ncol = db)
id <- rep(1:n,each=ni)
sigma <- 0.1*diag(ni)
omega <- 0.1*diag(db)


K<-1000
# Several Chains for the same iteration
M <- 1

alpha1 <- 0.7
alpha2 <- 0.4
seed0=44444
ylim <- c(0.3)

nsim <- 5
G<-2
col.names <- c("iteration", paste0("mu",1:G))
theta<-list(mu=mu)
theta0<-list(mu=mu0)
##  Simulation
x <- matrix(0,nrow=n*ni,ncol=nsim)
for (j in (1:nsim))
{
  seed <- j*seed0
  set.seed(seed)
  xj<-mixt.simulate(n,ni,A,B,mu,sigma,omega)
  x[,j] <- xj
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
	  em[,2:3] <- em[,2:3]+dem[dem$rep==j,2:3]
	}
}
em[,2:3] <- 1/nsim*em[,2:3]
em[,4]<-NULL
em_scaled <- em[rep(seq_len(nrow(em)), each=n),]
em_scaled$iteration = 0:(n*(K+1)-1)


## IEM
print('IEM')
nbr<-1
diem <- NULL
df.iem <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.iem_seq(x[,j], theta0, K,A,B,sigma,omega,id,nbr)
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df
}
# graphConvMC_new(diem, title="IEM 1R")
# diem[,2] <- diem[,2]^2
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


iem1$algo <- 'IEM'
em_scaled$algo <- 'EM'
variance <- NULL
variance <- rbind(em_scaled[0:K,],iem1)
colnames(variance) <- c("iteration","beta1","beta2","algo")
# graphConvMC2_new(variance, title="IEMs alpha=0.33",legend=TRUE)
plot_new(variance,legend=FALSE)



# ## IEM
# print('IEM')
# nbr<-50
# KR <- KNR*n/nbr
# diem <- NULL
# df.iem <- vector("list", length=nsim)
# for (j in (1:nsim))
# {
#   print(j)
#   df <- mixt.iem_seq(x[,j], theta0, KR, alph,nbr)
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

