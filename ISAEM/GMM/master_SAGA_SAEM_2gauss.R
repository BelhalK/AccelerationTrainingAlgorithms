source("mixtureAlgos.R")
source("mixtureFunctions.R")
theme_set(theme_bw())

#############################################
#####LS and precision plots (Master)####
#############################################


n <- 100
weight<-c(0.7, 0.3) 
mu<-c(0,4)
sigma<-c(1,1)*1


weight0<-c(.5,.5)
mu0<-c(1,2)
sigma0<-c(.5,2)
nb_r <- 10
KNR <- 500
K1 <-10
K <- 5000

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
  df <- mixt.ident(df)
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}
graphConvMC(dem, title="EM")

##  SAEM1 replacement vs EM
print('SAEM 10R')
nb_r <- 10
KR <- KNR*n/nb_r
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
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,table1r[i,2:7])
}


for (l in (0:(KR-1)))
{
  table1r[(l*nb_r+2):((l+1)*nb_r+1),2:7] <- Lr[l+1,]
}
table1r$iteration <- 1:(KR*nb_r+1)

##  SAEM1 replacement vs EM
print('SAEM 25R')
nb_r <- 25
KR <- KNR*n/nb_r
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
table2r$algo <- 'R'
table2r <- subset(table2r, select=c(7,1:8))
table2r[,8]<-NULL
table2r[,2:7] <- 1/nsim*table2r[,2:7]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,table2r[i,2:7])
}


for (l in (0:(KR-1)))
{
  table2r[(l*nb_r+2):((l+1)*nb_r+1),2:7] <- Lr[l+1,]
}
table2r$iteration <- 1:(KR*nb_r+1)


##  SAEM1 replacement vs EM
print('SAEM 50R')
nb_r <- 50
KR <- KNR*n/nb_r
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
table3r <- NULL
table3r <- diffr[diffr$sim==1,2:7]
for (j in (2:nsim))
{
  table3r <- table3r+diffr[diffr$sim==j,2:7]
}
table3r$iteration <- 0:KR
table3r$algo <- 'R'
table3r <- subset(table3r, select=c(7,1:8))
table3r[,8]<-NULL
table3r[,2:7] <- 1/nsim*table3r[,2:7]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,table3r[i,2:7])
}


for (l in (0:(KR-1)))
{
  table3r[(l*nb_r+2):((l+1)*nb_r+1),2:7] <- Lr[l+1,]
}
table3r$iteration <- 1:(KR*nb_r+1)

##  SAEM1 replacement vs EM
print('SAEM 75R')
nb_r <- 75
KR <- round(KNR*n/nb_r)-1
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
table4r <- NULL
table4r <- diffr[diffr$sim==1,2:7]
for (j in (2:nsim))
{
  table4r <- table4r+diffr[diffr$sim==j,2:7]
}
table4r$iteration <- 0:KR
table4r$algo <- 'R'
table4r <- subset(table4r, select=c(7,1:8))
table4r[,8]<-NULL
table4r[,2:7] <- 1/nsim*table4r[,2:7]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,table4r[i,2:7])
}


for (l in (0:(KR-1)))
{
  table4r[(l*nb_r+2):((l+1)*nb_r+1),2:7] <- Lr[l+1,]
}

table4r[49952:50001,2:7] <- Lr[666,]
table4r$iteration <- 1:50001

##  SAEM1 replacement vs EM
print('SAEM 65R')
nb_r <- 65
KR <- round(KNR*n/nb_r)-1
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
table6r <- NULL
table6r <- diffr[diffr$sim==1,2:7]
for (j in (2:nsim))
{
  table6r <- table6r+diffr[diffr$sim==j,2:7]
}
table6r$iteration <- 0:KR
table6r$algo <- 'R'
table6r <- subset(table6r, select=c(7,1:8))
table6r[,8]<-NULL
table6r[,2:7] <- 1/nsim*table6r[,2:7]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,table6r[i,2:7])
}


for (l in (0:(KR-1)))
{
  table6r[(l*nb_r+2):((l+1)*nb_r+1),2:7] <- Lr[l+1,]
}

table6r[49920:50001,2:7] <- Lr[768,]
table6r$iteration <- 1:50001


##  SAEM1 vs EM
print('SAEM NR')
diffnr <- NULL
for (j in (1:nsim))
{ 
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1(x[,j], theta0, KNR, K1, alpha=0.6, M)
  df <- mixt.ident(df)
  df <- df - df.em[[j]][1:(KNR+1),]
  df$iteration <- 0:KNR
  df$sim <- j
  diffnr <- rbind(diffnr,df)
}
diffnr[,2:7] <- diffnr[,2:7]^2
tablenr <- NULL
tablenr <- diffnr[diffnr$sim==1,2:7]
for (j in (2:nsim))
{
  tablenr <- tablenr+diffnr[diffnr$sim==j,2:7]
}
tablenr$iteration <- 0:KNR
tablenr$algo <- 'NR'
tablenr <- subset(tablenr, select=c(7,1:8))
tablenr[,8]<-NULL
tablenr[,2:7] <- 1/nsim*tablenr[,2:7]

Lnr <- NULL
for (i in (2:(KNR+1)))
{
  Lnr <- rbind(Lnr,tablenr[i,2:7])
}

for (ind in (0:(KNR-1)))
{
  tablenr[(ind*n+2):((ind+1)*n+1),2:7] <- Lnr[ind+1,]
}
tablenr$iteration <- 1:(KNR*n+1)

# tablenr[1,2:7] <- tablenr[2,2:7]
table1r$algo <- '10R'
table2r$algo <- '25R'
table3r$algo <- '50R'
table4r$algo <- '75R'
table6r$algo <- '65R'
tablenr$algo <- 'NR'

variance <- rbind(table1r[1:5000,],table2r[1:5000,],table3r[1:5000,],table4r[1:5000,],table6r[1:5000,],tablenr[1:5000,]) #10replacement
variance <- rbind(table6r[1:5000,],table3r[1:5000,],tablenr[1:5000,]) #10replacement
# var <- graphConvMC2(variance, title="ALGO - EM (same complexity)",legend=TRUE)
# ggsave('conv_100sim.png',var)
graphConvMC2(variance, title="ALGO - EM (same complexity)",legend=TRUE)

precision <- NULL
precision <- variance[1:6,2:7]

# precision[1,] <- 1/10000*colSums(tail(table1r[,2:7],10000))
# precision[2,] <- 1/10000*colSums(tail(table2r[,2:7],10000))
# precision[3,] <- 1/10000*colSums(tail(table3r[,2:7],10000))
# precision[4,] <- 1/10000*colSums(tail(table6r[,2:7],10000))
# precision[5,] <- 1/10000*colSums(tail(table4r[,2:7],10000))
# precision[6,] <- 1/10000*colSums(tail(tablenr[,2:7],10000))
precision[1,] <- table1r[40000,2:7]
precision[2,] <- table2r[40000,2:7]
precision[3,] <- table3r[40000,2:7]
precision[4,] <- table6r[40000,2:7]
precision[5,] <- table4r[40000,2:7]
precision[6,] <- tablenr[40000,2:7]
precision$iteration <- 0:5
precision <- subset(precision, select=c(7,1:7))
precision[,8]<-NULL
precision$algo <- 'SAEM'

prec <- graphConvMC2(precision, title="precision = f(%R)",legend=TRUE)
ggsave('precision3.png',prec)

