source("mixtureAlgos.R")
source("mixtureFunctions.R")
theme_set(theme_bw())

n <- 100
weight<-c(0.4, 0.3,0.3) 
mu<-c(0,2,4)
sigma<-c(1,1,1)*0.6


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

nsim <- 5
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
dem <- NULL
df.em <- vector("list", length=nsim)
for (j in (1:nsim))
{
  df <- mixt.em(x[,j], theta, K)
  df <- mixt.ident3(df)
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}
graphConvMC_new(dem, title="EM")


##  SAEM2
diff <- NULL
for (j in (1:nsim))
{
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1(x[,j], theta0, K, K1, M=1, alpha=0.6)
  df <- mixt.ident3(df)
  df <- df - df.em[[j]]
  df$iteration <- 0:K
  df$rep <- j
  diff <- rbind(diff,df)
}
graphConvMC_new(diff, title="SAEM - EM", ylim=ylim)




##  SAEMRP
print('SAEM 10R')
nb_r <- 10
KR <- KNR*n/nb_r
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M=1, nb_r)
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
table1r$algo <- '10R'
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
# graphConvMC_new(diffr, title="1OR")
graphConvMC2_new(table1r, title="10RALGO - EM (same complexity)",legend=TRUE)

##  SAEMRP
print('SAEM 50R')
nb_r <- 50
KR <- KNR*n/nb_r
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M=1, nb_r)
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
table2r$algo <- '50R'
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
# graphConvMC_new(diffr, title="5OR")
graphConvMC2_new(table2r, title="50RALGO - EM (same complexity)",legend=TRUE)

##  SAEM75R
print('SAEM 75R')
nb_r <- 75
KR <- KNR*n/nb_r
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M=1, nb_r)
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
table3r$algo <- '75R'
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
table3r[49950:50001,2:10] <- Lr[665,]
table3r$iteration <- 1:(KR*nb_r+1)
# graphConvMC_new(diffr, title="5OR")
graphConvMC2_new(table3r, title="75RALGO - EM (same complexity)",legend=TRUE)


##  SAEMNR
print('SAEM NR')

KR <- KNR
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1(x[,j], theta0, KR, K1, alpha=0.6, M=1)
  df <- mixt.ident3(df)
  df <- df - df.em[[j]][1:(KR+1),]
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
tablenr$iteration <- 0:KR
tablenr$algo <- 'NR'
tablenr <- subset(tablenr, select=c(10,1:11))
tablenr[,11]<-NULL
tablenr[,2:10] <- 1/nsim*tablenr[,2:10]
Lr <- NULL
for (i in (2:(KR+1)))
{
  Lr <- rbind(Lr,tablenr[i,2:10])
}


for (l in (0:(KR-1)))
{
  tablenr[(l*n+2):((l+1)*n+1),2:10] <- Lr[l+1,]
}
tablenr$iteration <- 1:(KR*n+1)

# graphConvMC_new(diffr, title="5OR")
graphConvMC2_new(tablenr, title="NRRALGO - EM (same complexity)",legend=TRUE)

table1r$algo <- '10R'
table2r$algo <- '50R'
table3r$algo <- '75R'
tablenr$algo <- 'NR'

variance <- NULL
variance <- rbind(table1r,table2r,table3r,tablenr) #10replacement
graphConvMC2_new(variance, title="ALGO - EM (same complexity)",legend=TRUE)
graphConvMC2_new(rbind(table2r, tablenr), title="50 vs NR",legend=TRUE)


# pdfem <- NULL
# df.pdfem <- vector("list", length=nsim)
# for (j in (1:nsim))
# {
#   param <- df.em[[j]]
#   data<- x[,j]
#   df <- mixt.pdff(param,data,K,n)
#   df$rep <- j
#   pdfem <- rbind(pdfem,df)
#   df$rep <- NULL
#   df.pdfem[[j]] <- df
# }
# graphConvMC_new(pdfem, title="EM")


# pdfsaem <- NULL
# df.pdfsaem <- vector("list", length=nsim)
# for (j in (1:nsim))
# {
#   param <- diff[diff$rep==j,1:10]
#   data<- x[,j]
#   df <- mixt.pdff(param,data,K,n)
#   df$rep <- j
#   pdfsaem <- rbind(pdfsaem,df)
#   df$rep <- NULL
#   df.pdfsaem[[j]] <- df
# }
# graphConvMC_new(pdfsaem, title="EM")

# pdfsaemr <- NULL
# df.pdfsaemr <- vector("list", length=nsim)
# for (j in (1:nsim))
# {
#   param <- diffr[diffr$rep==j,1:10]
#   data<- x[,j]
#   df <- mixt.pdff(param,data,K,n)
#   df$rep <- j
#   pdfsaemr <- rbind(pdfsaemr,df)
#   df$rep <- NULL
#   df.pdfsaemr[[j]] <- df
# }
# graphConvMC_new(pdfsaemr, title="EM")

