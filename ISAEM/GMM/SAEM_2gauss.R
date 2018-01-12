source("mixtureAlgos.R")
source("mixtureFunctions.R")
theme_set(theme_bw())


n <- 100
weight<-c(0.7, 0.3) 
mu<-c(0,4)
sigma<-c(1,1)*1


weight0<-c(.5,.5)
mu0<-c(1,2)
sigma0<-c(.5,2)
nb_r <- 10
K1 <-50
K <- 100

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
  df <- mixt.ident(df)
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}
graphConvMC(dem, title="EM")

pdf <- NULL
param <- df.em[[1]]
data<- x[,1]
pdf <- mixt.pdf(param,data,K,n)
pdftest <- mixt.pdftest(param,data,K,n)
plot(pdf)
plot(pdftest)

##  SAEM2
diff <- NULL
for (j in (1:nsim))
{
  seed <- j*seed0
  set.seed(seed)
  df <- mixt.saem1(x[,j], theta0, K, K1, M=1, alpha=0.6)
  df <- mixt.ident(df)
  df <- df - df.em[[j]]
  df$iteration <- 0:K
  df$rep <- j
  diff <- rbind(diff,df)
}
graphConvMC(diff, title="SAEM - EM", ylim=ylim)
