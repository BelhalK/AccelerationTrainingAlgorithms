require(ggplot2)
require(gridExtra)
require(reshape2)
library(rlist)

source("utils/algos.R")
source("utils/func.R")
source("utils/plots.R")
theme_set(theme_bw())
options(digits = 2)

n <- 1000
nb.epochs <- 15
K <- n*nb.epochs
nsim=10

weight<-c(0.2, 0.8)
mean <- 0.5
mu<-c(mean,-mean)
sigma<-c(1,1)*1


weight0<-weight
mean0 <- 1.1
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


nbr<-1
nbr10 <- n/200
# nbr10 <- n/10
nbr50 <- n/2

print('EM')
dem <- NULL
df.em <- vector("list", length=nsim)

dsaem <- NULL
df.saem <- vector("list", length=nsim)
Kem <- nbr*K/n

diemseq <- NULL
df.iemseq <- vector("list", length=nsim)

disaem <- NULL
df.isaem <- vector("list", length=nsim)
disaem10 <- NULL
df.isaem10 <- vector("list", length=nsim)
disaem50 <- NULL
df.isaem50 <- vector("list", length=nsim)

kiter = 1:K


nb.chains <- 30

for (j in (1:nsim))
{

  print(j)
  seed <- j*seed0
  set.seed(seed)
  ML <- mls[[j]]
  print("ML calculation done")

  df <- mixt.em(x[,j], theta0, nb.epochs)
  df[,2:7] <- (df[,2:7] - ML[1:(Kem+1),2:7])^2
  # df[,2:7] <- (df[,2:7] - ML[,2:7])^2
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
  print('em done')

  df <- mixt.saem(x[,j],theta0, nb.epochs, K1=Kem/2, alpha=0.6, M=1)
  df[,2:7] <- (df[,2:7] - ML[,2:7])^2
  df$rep <- j
  dsaem <- rbind(dsaem,df)
  df$rep <- NULL
  df.saem[[j]] <- df
  print('saem done')


  df <- mixt.iem.seq(x[,j], theta0, nb.epochs*n/nbr,nbr)
  df[,2:7] <- (df[,2:7] - ML[,2:7])^2
  df$rep <- j
  diemseq <- rbind(diemseq,df)
  df$rep <- NULL
  df.iemseq[[j]] <- df
  print('iemseq done')

  df <- mixt.isaem(x[,j],theta0, nb.epochs*n/nbr, K1=K/2, alpha=0.6, M=nb.chains,nbr)
  df[,2:7] <- (df[,2:7] - ML[,2:7])^2
  df$rep <- j
  disaem <- rbind(disaem,df)
  df$rep <- NULL
  df.isaem[[j]] <- df
  print('isaem done')

  df <- mixt.isaem(x[,j],theta0, nb.epochs*n/nbr10, K1=K/2, alpha=0.6, M=nb.chains,nbr10)
  df[,2:7] <- (df[,2:7] - ML[,2:7])^2
  df$rep <- j
  disaem10 <- rbind(disaem10,df)
  df$rep <- NULL
  df.isaem10[[j]] <- df
  print('isaem10 done')

  df <- mixt.isaem(x[,j],theta0, nb.epochs*n/nbr50, K1=K/2, alpha=0.6, M=nb.chains,nbr50)
  df[,2:7] <- (df[,2:7] - ML[,2:7])^2
  df$rep <- j
  disaem50 <- rbind(disaem50,df)
  df$rep <- NULL
  df.isaem50[[j]] <- df
  print('isaem50 done')



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

saem <- NULL
saem <- dsaem[dsaem$rep==1,]

if (nsim>2) {
   for (j in (2:nsim))
  {
    saem[,2:7] <- saem[,2:7]+dsaem[dsaem$rep==j,2:7]
  }
}
saem[,2:7] <- 1/nsim*saem[,2:7]
saem[,9]<-NULL



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



isaem <- NULL
isaem <- disaem[disaem$rep==1,]


if (nsim>2) {
    for (j in (2:nsim))
  {
    isaem[,2:7] <- isaem[,2:7]+disaem[disaem$rep==j,2:7]
  }
}

isaem[,2:7] <- 1/nsim*isaem[,2:7]
isaem[,9]<-NULL


isaem10 <- NULL
isaem10 <- disaem10[disaem10$rep==1,]


if (nsim>2) {
    for (j in (2:nsim))
  {
    isaem10[,2:7] <- isaem10[,2:7]+disaem10[disaem10$rep==j,2:7]
  }
}

isaem10[,2:7] <- 1/nsim*isaem10[,2:7]
isaem10[,9]<-NULL

isaem50 <- NULL
isaem50 <- disaem50[disaem50$rep==1,]


if (nsim>2) {
    for (j in (2:nsim))
  {
    isaem50[,2:7] <- isaem50[,2:7]+disaem50[disaem50$rep==j,2:7]
  }
}

isaem50[,2:7] <- 1/nsim*isaem50[,2:7]
isaem50[,9]<-NULL

em$algo <- 'EM'
iemseq$algo <- 'IEM seq'
saem$algo <- 'saem'
isaem$algo <- 'isaem'
isaem10$algo <- 'isaem10'
isaem50$algo <- 'isaem50'

em$rep <- NULL
iemseq$rep <- NULL
saem$rep <- NULL
isaem$rep <- NULL
isaem10$rep <- NULL
isaem50$rep <- NULL


iemseq$iteration <- iemseq$iteration*nbr/n
iemseq$algo <- 'IEM'
isaem$iteration <- isaem$iteration*nbr/n
isaem10$iteration <- isaem10$iteration*nbr10/n
isaem50$iteration <- isaem50$iteration*nbr50/n


variance <- rbind(isaem[,c(1,4,8)],
                  isaem10[,c(1,4,8)],
                  isaem50[,c(1,4,8)],
                  iemseq[,c(1,4,8)],
                  em[,c(1,4,8)],
                  saem[,c(1,4,8)])

graphConvMC2_new(variance, title="",legend=TRUE)

save.image("gmm.RData")
write.csv(variance, file = "notebooks/singlerun.csv")



# ### PER EPOCH
# epochs = seq(1, K, by=n/nbr)
# em_ep <- em[1:(nbr*K/n),]
# em_ep$iteration <- 1:(nbr*K/n)
# saem_ep <- saem[1:(nbr*K/n),]
# saem_ep$iteration <- 1:(nbr*K/n)
# iem_ep <- iemseq[epochs,]
# iem_ep$iteration <- 1:(K/n)
# isaem_ep <- isaem[epochs,]
# isaem_ep$iteration <- 1:(K/n)


# start =1000
# end = K


# em_ep <- em
# em_ep$iteration <- n*em$iteration

# saem_ep <- saem
# saem_ep$iteration <- n*saem$iteration


# # variance <- rbind(iemseq[start:end,c(1,4,8)],
# #                   em_ep[2:length(epochs),c(1,4,8)],
# #                   saem_ep[2:length(epochs),c(1,4,8)],
# #                   isaem[start:end,c(1,4,8)])

# # graphConvMC2_new(variance, title="",legend=TRUE)



# start =1000
# end = K

# testiemseq <- iemseq
# testisaem <- isaem
# testem_ep <- em_ep
# testsaem_ep <- saem_ep


# testiemseq$iteration <- testiemseq$iteration/n
# testiemseq$algo <- 'IEM'
# testisaem$iteration <- testisaem$iteration/n
# testem_ep$iteration <- testem_ep$iteration/n
# testsaem_ep$iteration <- testsaem_ep$iteration/n



# variance <- rbind(testisaem[start:end,c(1,4,8)],
#                   testiemseq[start:end,c(1,4,8)],
#                   testem_ep[2:length(epochs),c(1,4,8)],testsaem_ep[2:length(epochs),c(1,4,8)])

# graphConvMC2_new(variance, title="",legend=TRUE)

