source("mixtureAlgos.R")
source("mixtureFunctions.R")
theme_set(theme_bw())

n <- 200
weight<-c(0.4, 0.3,0.3) 
mu<-c(0,2,4)
sigma<-c(0.3,0.3,0.3)

weight0<-c(0.1, 0.3,0.6) 
mu0<-c(1,2,3)
sigma0<-c(0.2,0.1,0.3)

KNR <- 500
K1 <-10
K <- 5000
# Several Chains for the same iteration
# M <- 10

alpha1 <- 0.7
alpha2 <- 0.4
seed0=44444

# ylim <- c(0.15, 0.5, 0.4)
ylim <- c(0.1, 0.3, 0.3)

nsim <- 30
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
ML <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.em(x[,j], theta0, K)
  # df <- mixt.ident3(df)
  df$sim <- j
  dem <- rbind(dem,df)
  df$sim <- NULL
  df.em[[j]] <- df
  ML[[j]] <- df
  ML[[j]][1:(K+1),]<- df[(K+1),]
}
graphConvMC_new(dem, title="EM")

##  SAEMRP
print('SAEM 10R')
nb_r <- n/10
KR <- K
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M=1, nb_r)
  df <- df - ML[[j]]
  df$iteration <- 0:KR
  df$sim <- j
  diffr <- rbind(diffr,df)
}
# graphConvMC_new(diffr, title="EM")
diffr[,2:10] <- diffr[,2:10]^2
table1r <- NULL
table1r <- diffr[diffr$sim==1,2:10]
for (j in (2:nsim))
{
  table1r <- table1r+diffr[diffr$sim==j,2:10]
}
table1r$iteration = seq(0, K, by=1)
table1r <- subset(table1r, select=c(10,1:9))
table1r[,2:10] <- 1/nsim*table1r[,2:10]


##  SAEMRP
print('SAEM 50R')
nb_r <- n/2
KR <- K
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M=1, nb_r)
  df <- df - ML[[j]]
  df$iteration <- 0:KR
  df$sim <- j
  diffr <- rbind(diffr,df)
}
diffr[,2:10] <- diffr[,2:10]^2
table2r <- NULL
table2r <- diffr[diffr$sim==1,2:10]
for (j in (2:nsim))
{
  table2r <- table2r+diffr[diffr$sim==j,2:10]
}
table2r$iteration = seq(0, 5*K, by=5)
table2r <- subset(table2r, select=c(10,1:9))
table2r[,2:10] <- 1/nsim*table2r[,2:10]
table2r <- table2r[rep(seq_len(nrow(table2r)), each=5),]


print('SAEM 100R')
nb_r <- n
KR <- K
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M=1, nb_r)
  df <- df - ML[[j]]
  df$iteration <- 0:KR
  df$sim <- j
  diffr <- rbind(diffr,df)
}

diffr[,2:10] <- diffr[,2:10]^2
tablenr <- NULL
tablenr <- diffr[diffr$sim==1,2:10]
for (j in (2:nsim))
{
  tablenr <- tablenr+diffr[diffr$sim==j,2:10]
}
tablenr$iteration = seq(0, 10*K, by=10)
tablenr <- subset(tablenr, select=c(10,1:9))
tablenr[,2:10] <- 1/nsim*tablenr[,2:10]
tablenr <- tablenr[rep(seq_len(nrow(tablenr)), each=10),]

# graphConvMC_new(diffr, title="5OR")
# graphConvMC2_new(tablenr, title="NRRALGO - EM (same complexity)",legend=TRUE)

# tablenr[1,2:7] <- tablenr[2,2:7]
table1r$algo <- '10R'
table2r$algo <- '50R'
tablenr$algo <- 'NR'

variance2 <- rbind(table1r[1:K,],table2r[1:K,],tablenr[1:K,]) 
graphConvMC2_new(variance2, title="GMM 2",legend=TRUE)
variance <- rbind(table1r[11:K,],table2r[11:K,],tablenr[11:K,]) 
graphConvMC2_new(variance, title="GMM2 zoom",legend=TRUE)
# var <- graphConvMC2(variance, title="ALGO - EM (same complexity)",legend=TRUE)
# ggsave('conv_100sim.png',var)


