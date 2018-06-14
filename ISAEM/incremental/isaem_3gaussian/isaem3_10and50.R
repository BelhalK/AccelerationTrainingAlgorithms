source("mixtureAlgos.R")
source("mixtureFunctions.R")
theme_set(theme_bw())

#############################################
#####LS and precision plots (Master)####
#############################################

n <- 100
weight<-c(0.5, 0.2,0.3) 
mu<-c(0,2,4)
sigma<-c(1,0.5,1)*1



K1 <-10
K <- 500
KNR <- 50

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
  df <- mixt.ident3(df)
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}
# graphConvMC(dem, title="EM")

##  SAEM1 replacement vs EM
print('SAEM original 50R')
nb_r <- 50
KR <- KNR*n/nb_r
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident3(df)
  df <- df - df.em[[j]][1:(KR+1),]
  df$iteration <- 0:KR
  df$sim <- j
  diffr <- rbind(diffr,df)
}
diffr[,2:10] <- diffr[,2:10]^2
table1r <- NULL
table1r <- diffr[diffr$sim==1,2:10]
for (j in (2:nsim))
{
  table1r <- table1r+diffr[diffr$sim==j,2:10]
}
table1r$iteration <- 0:KR
table1r$algo <- 'R'
table1r <- subset(table1r, select=c(10,1:11))
table1r[,11]<-NULL
table1r[,2:10] <- 1/nsim*table1r[,2:10]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,table1r[i,2:10])
}


for (l in (0:(KR-1)))
{
  table1r[(l*nb_r+2):((l+1)*nb_r+1),2:10] <- Lr[l+1,]
}
table1r$iteration <- 1:(KR*nb_r+1)




##  SAEM1new replacement vs EM
print('SAEM NEW 50R')
nb_r <- 50
KR <- KNR*n/nb_r
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.isaem3(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident3(df)
  df <- df - df.em[[j]][1:(KR+1),]
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
table2r$iteration <- 0:KR
table2r$algo <- 'R'
table2r <- subset(table2r, select=c(10,1:11))
table2r[,11]<-NULL
table2r[,2:10] <- 1/nsim*table2r[,2:10]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,table2r[i,2:10])
}


for (l in (0:(KR-1)))
{
  table2r[(l*nb_r+2):((l+1)*nb_r+1),2:10] <- Lr[l+1,]
}
table2r$iteration <- 1:(KR*nb_r+1)

##  SAEM1new replacement vs EM
print('SAEM NEW 50R')
nb_r <- 10
KR <- KNR*n/nb_r
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident3(df)
  df <- df - df.em[[j]][1:(KR+1),]
  df$iteration <- 0:KR
  df$sim <- j
  diffr <- rbind(diffr,df)
}
diffr[,2:10] <- diffr[,2:10]^2
table3r <- NULL
table3r <- diffr[diffr$sim==1,2:10]
for (j in (2:nsim))
{
  table3r <- table3r+diffr[diffr$sim==j,2:10]
}
table3r$iteration <- 0:KR
table3r$algo <- 'R'
table3r <- subset(table3r, select=c(10,1:11))
table3r[,11]<-NULL
table3r[,2:10] <- 1/nsim*table3r[,2:10]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,table3r[i,2:10])
}


for (l in (0:(KR-1)))
{
  table3r[(l*nb_r+2):((l+1)*nb_r+1),2:10] <- Lr[l+1,]
}
table3r$iteration <- 1:(KR*nb_r+1)

##  SAEM1new replacement vs EM
print('SAEM NEW 50R')
nb_r <- 10
KR <- KNR*n/nb_r
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.isaem3(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident3(df)
  df <- df - df.em[[j]][1:(KR+1),]
  df$iteration <- 0:KR
  df$sim <- j
  diffr <- rbind(diffr,df)
}
diffr[,2:10] <- diffr[,2:10]^2
table4r <- NULL
table4r <- diffr[diffr$sim==1,2:10]
for (j in (2:nsim))
{
  table4r <- table4r+diffr[diffr$sim==j,2:10]
}
table4r$iteration <- 0:KR
table4r$algo <- 'R'
table4r <- subset(table4r, select=c(10,1:11))
table4r[,11]<-NULL
table4r[,2:10] <- 1/nsim*table4r[,2:10]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,table4r[i,2:10])
}


for (l in (0:(KR-1)))
{
  table4r[(l*nb_r+2):((l+1)*nb_r+1),2:10] <- Lr[l+1,]
}
table4r$iteration <- 1:(KR*nb_r+1)


##  SAEM1 vs EM
print('SAEM NR')
diffnr <- NULL
for (j in (1:nsim))
{ 
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1(x[,j], theta0, KNR, K1, alpha=0.6, M)
  df <- mixt.ident3(df)
  df <- df - df.em[[j]][1:(KNR+1),]
  df$iteration <- 0:KNR
  df$sim <- j
  diffnr <- rbind(diffnr,df)
}
diffnr[,2:10] <- diffnr[,2:10]^2
tablenr <- NULL
tablenr <- diffnr[diffnr$sim==1,2:10]
for (j in (2:nsim))
{
  tablenr <- tablenr+diffnr[diffnr$sim==j,2:10]
}
tablenr$iteration <- 0:KNR
tablenr$algo <- 'NR'
tablenr <- subset(tablenr, select=c(10,1:11))
tablenr[,11]<-NULL
tablenr[,2:10] <- 1/nsim*tablenr[,2:10]

Lnr <- NULL
for (i in (2:(KNR+1)))
{
  Lnr <- rbind(Lnr,tablenr[i,2:10])
}


for (ind in (0:(KNR-1)))
{
  tablenr[(ind*n+2):((ind+1)*n+1),2:10] <- Lnr[ind+1,]
}
tablenr$iteration <- 1:(KNR*n+1)

table1r$algo <- '50R'
table2r$algo <- '50R 3'
table3r$algo <- '10R'
table4r$algo <- '10R 3'
tablenr$algo <- 'NR'

variance <- NULL
# variance <- rbind(table1r,tablenr)
variance <- rbind(table1r,table2r,table3r,table4r,tablenr)
# variance <- rbind(table1r,table2r,table3r,table4r,tablenr)
var <- graphConvMC2_new(variance, title="ALGO - EM (same complexity)",legend=TRUE)
ggsave('isaem3_1050.png',var)

