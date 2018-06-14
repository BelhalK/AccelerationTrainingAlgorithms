source("mixtureAlgos.R")
source("mixtureFunctions.R")
theme_set(theme_bw())

#############################################
##### ISAEM with Bouchard techniques to sample the individual
#############################################

n <- 100
weight<-c(0.3, 0.7) 
mu<-c(0,2)
sigma<-c(1,0.3)*1



K1 <-10
K <- 500
KNR <- 250

alpha1 <- 0.7
alpha2 <- 0.4
seed0=44444


# ylim <- c(0.15, 0.5, 0.4)
ylim <- c(0.1, 0.3, 0.3)

M <- 1
nsim <- 50
nb_r <- 50
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
# graphConvMC_new(dem, title="EM")



##  SAEM1 replacement vs EM
print('SAEM original 50R')

KR<-KNR
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1_replace1(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident(df)
  for (k in (0:KR))
  {
    df[k+1,2:7] <- (df[k+1,2:7] - df.em[[j]][(K),2:7])^2
  }
  # df <- df - df.em[[j]][1:(KR+1),]
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
table1r$algo <- 'R'
table1r <- subset(table1r, select=c(7,1:8))
table1r[,8]<-NULL
table1r[,2:7] <- 1/nsim*table1r[,2:7]



#PISAEM
print('SAEM NEW 50R')

KR<-KNR
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.isaem2(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident(df)
  for (k in (0:KR))
  {
    df[k+1,2:7] <- (df[k+1,2:7] - df.em[[j]][(K),2:7])^2
  }
  # df <- df - df.em[[j]][1:(KR+1),]
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

#ISAEM with grades
print('SAEM NEW 50R')

KR<-KNR
diffr <- NULL

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.isaem3(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
  df <- mixt.ident(df)
  for (k in (0:KR))
  {
    df[k+1,2:7] <- (df[k+1,2:7] - df.em[[j]][(K),2:7])^2
  }
  # df <- df - df.em[[j]][1:(KR+1),]
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
  for (k in (0:KNR))
  {
    df[k+1,2:7] <- (df[k+1,2:7] - df.em[[j]][(K),2:7])^2
  }
  # df <- df - df.em[[j]][1:(KNR+1),]
  df$iteration <- 0:KNR
  df$sim <- j
  diffnr <- rbind(diffnr,df)
}

# names(diffnr)[8] <- paste("rep")
# graphConvMC_new(diffnr, title="ALGO - EM (same complexity)")
# names(diffnr)[11] <- paste("sim")

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

tablenr <- tablenr[rep(1:nrow(tablenr),each=n/nb_r),]
tablenr <- tablenr[1:nrow(table1r),]

tablenr$iteration <- table1r$iteration




# ##  SAEM1new replacement vs EM
# print('SAEM NEW 50R')

# KR<-KNR
# diffr <- NULL

# for (j in (1:nsim))
# {
#   print(j)
#   seed <- j*seed0
#   set.seed(seed)
#   df <- mixt.isaem_bouchard(x[,j], theta0, KR, K1, alpha=0.6, M, nb_r)
#   df <- mixt.ident3(df)
#   for (k in (0:KR))
#   {
#     df[k+1,2:7] <- (df[k+1,2:7] - df.em[[j]][(K),2:7])^2
#   }
#   # df <- df - df.em[[j]][1:(KR+1),]
#   df$iteration <- 0:KR
#   df$sim <- j
#   diffr <- rbind(diffr,df)
# }

# # graphConvMC_new(diffr, title="EM")

# diffr[,2:7] <- diffr[,2:7]^2
# table2r <- NULL
# table2r <- diffr[diffr$sim==1,2:7]
# for (j in (2:nsim))
# {
#   table2r <- table2r+diffr[diffr$sim==j,2:7]
# }
# table2r$iteration <- 0:KR
# table2r$algo <- 'R'
# table2r <- subset(table2r, select=c(7,1:8))
# table2r[,8]<-NULL
# table2r[,2:7] <- 1/nsim*table2r[,2:7]



table1r$algo <- 'ISAEM rand'
# table2r$algo <- 'Bouch'
table3r$algo <- 'PISAEM'
table4r$algo <- 'GRADES'
tablenr$algo <- 'NR'

variance <- NULL
# variance <- rbind(table2r)
# variance <- rbind(table3r,table4r)
# variance <- rbind(table1r,table3r,tablenr)
variance <- rbind(table1r,table3r,table4r,tablenr)
# variance <- rbind(table1r,table2r,table3r,tablenr)
# variance <- rbind(table4r,table2r,table3r)
# variance <- rbind(table1r,table2r,table3r,table4r,tablenr)

# variance <- rbind(table4r,table2r,table3r, tablenr)


var <- melt(variance, id.var = c('iteration','algo'), na.rm = TRUE)
var <- graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)

