source("mixtureAlgos.R")
source("mixtureFunctions.R")
theme_set(theme_bw())

#############################################
#####LS and precision plots (Master)####
#############################################


n <- 500
weight<-c(0.7, 0.3) 
mu<-c(0,1)
sigma<-c(0.2,0.2)*1


# weight0<-c(.5,.5)
# mu0<-c(1,2)
# sigma0<-c(.5,2)

KNR <- 500
K1 <-10
K <- 5000

alpha1 <- 0.7
alpha2 <- 0.4
seed0=44444


# ylim <- c(0.15, 0.5, 0.4)
ylim <- c(0.1, 0.3, 0.3)

M <- 1
nsim <- 15
#
G<-length(mu)
col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
theta<-list(p=weight,mu=mu,sigma=sigma)
# theta0<-list(p=weight0,mu=mu0,sigma=sigma0)
theta0<-theta


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
  df <- mixt.em(x[,j], theta, K)
  df <- mixt.ident(df)
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}
# graphConvMC(dem, title="EM")

##  SAEM1 replacement vs EM
print('SAEM 10R')
nb_r <- n/10
KR <- 100
diffr <- NULL
for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident(df)
  df <- df - df.em[[j]][1:(KR+1),]
  df$iteration <- 0:KR
  df$sim <- j
  diffr <- rbind(diffr,df)
}
diffr[,2:7] <- diffr[,2:7]^2
table1r <- NULL
table1r <- diffr[diffr$sim==1,2:7]
for (j in (2:nsim))
{
  table1r <- table1r+diffr[diffr$sim==j,2:7]
}
table1r$iteration <- 0:KR
table1r$algo <- '10R'
table1r <- subset(table1r, select=c(7,1:8))
table1r[,8]<-NULL
table1r[,2:7] <- 1/nsim*table1r[,2:7]

print('SAEM 100R')
nb_r <- n
KR <- 100
diffr <- NULL
for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident(df)
  df <- df - df.em[[j]][1:(KR+1),]
  df$iteration <- 0:KR
  df$sim <- j
  diffr <- rbind(diffr,df)
}
diffr[,2:7] <- diffr[,2:7]^2
tablenr <- NULL
tablenr <- diffr[diffr$sim==1,2:7]
for (j in (2:nsim))
{
  tablenr <- tablenr+diffr[diffr$sim==j,2:7]
}
tablenr$iteration <- 0:KR
tablenr$algo <- 'NR'
tablenr <- subset(tablenr, select=c(7,1:8))
tablenr[,8]<-NULL
tablenr[,2:7] <- 1/nsim*tablenr[,2:7]

tablenr_scaled <- tablenr[rep(seq_len(nrow(tablenr)), each=100/10),]
tablenr_scaled$iteration = 1:(10*(KR+1))


print('SAEM 100R')
nb_r <- n/2
KR <- 100
diffr <- NULL
for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident(df)
  df <- df - df.em[[j]][1:(KR+1),]
  df$iteration <- 0:KR
  df$sim <- j
  diffr <- rbind(diffr,df)
}
diffr[,2:7] <- diffr[,2:7]^2
table2r <- NULL
table2r <- diffr[diffr$sim==1,2:7]
for (j in (2:nsim))
{
  table2r <- table2r+diffr[diffr$sim==j,2:7]
}
table2r$iteration <- 0:KR
table2r$algo <- 'NR'
table2r <- subset(table2r, select=c(7,1:8))
table2r[,8]<-NULL
table2r[,2:7] <- 1/nsim*table2r[,2:7]


table2r_scaled <- table2r[rep(seq_len(nrow(table2r)), each=5),]
tablenr_scaled <- tablenr[rep(seq_len(nrow(tablenr)), each=10),]
table2r_scaled$iteration = 1:(5*(KR+1))
tablenr_scaled$iteration = 1:(10*(KR+1))



# tablenr[1,2:7] <- tablenr[2,2:7]
table1r$algo <- '10R'
table2r$algo <- '50R'
tablenr_scaled$algo <- 'NR'

variance <- rbind(table1r[1:KR,],table2r[1:KR,],tablenr[1:KR,]) #10replacement
# var <- graphConvMC2(variance, title="ALGO - EM (same complexity)",legend=TRUE)
# ggsave('conv_100sim.png',var)
graphConvMC2(variance, title="ALGO - EM (same complexity)",legend=TRUE)
