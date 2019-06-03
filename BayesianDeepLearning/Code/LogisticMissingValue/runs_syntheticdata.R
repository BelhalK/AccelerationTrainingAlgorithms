library(MASS)
library(ggplot2)
library(reshape2)

source("R/miss.saem.R")
source("R/louis_lr_saem.R")
source("R/log_reg.R")
source("R/likelihood_saem.R")

#Generate dataset
N <- 100  # number of subjects
p <- 3     # number of explanatory variables
mu.star <- rep(0,p)  # mean of the explanatory variables
Sigma.star <- diag(rep(1,p)) # covariance
beta.star <- c(1, 1,  0) # coefficients
beta0.star <- 0 # intercept
beta.true = c(beta0.star,beta.star)
X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) +
              matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
y <- as.numeric(runif(N)<p1)

# Generate missingness
p.miss <- 0.10
patterns <- runif(N*p)<p.miss #missing completely at random
X.obs <- X.complete
X.obs[patterns] <- NA


batchsize = 1
nb.epochs <-10
nb.iter <- N/batchsize*nb.epochs

# SAEM
list.saem = miss.saem(X.obs,y,maxruns=nb.iter,k1=0,ll_obs_cal=FALSE, algo = "saem")
#MCEM
list.mcem = miss.saem(X.obs,y,maxruns=nb.iter,ll_obs_cal=FALSE, algo = "mcem")
#Incremental MCEM
list.imcem = miss.saem(X.obs,y,maxruns=nb.iter,ll_obs_cal=FALSE, algo = "imcem", batchsize= batchsize)



#PLOTS
dim = 1
imcem = list.imcem$mu[dim,]
saem = list.saem$mu[dim,]
mcem = list.mcem$mu[dim,]

imcem = list.imcem$seqbeta[dim,]
saem = list.saem$seqbeta[dim,]
mcem = list.mcem$seqbeta[dim,]

x = 1:length(imcem)
df <- data.frame(x,saem,mcem, imcem)
ggplot(data=df)+
  geom_line(mapping=aes(y=saem,x= x,color="saem"),size=0.5 ) +
  geom_line(mapping=aes(y=mcem,x= x,color="mcem"),size=0.5) +
  geom_line(mapping=aes(y=imcem,x= x,color="imcem"),size=0.5) +
  scale_color_manual(values = c(
    'saem' = 'darkblue','mcem' = 'red','imcem' = 'black')) +
  labs(color = 'Algo')+ ylab("beta")


#PER EPOCHS
epochs = seq(1,nb.iter,N/batchsize)
x = 2:(nb.epochs+1)
saem.ep <- saem[x]
mcem.ep <- mcem[x]
imcem.ep <- imcem[(epochs+1)]
df <- data.frame(x,saem.ep,mcem.ep,imcem.ep)

ggplot(data=df)+
  geom_line(mapping=aes(y=saem.ep,x= x,color="saem"),size=0.5 ) +
  geom_line(mapping=aes(y=mcem.ep,x= x,color="mcem"),size=0.5) +
  geom_line(mapping=aes(y=imcem.ep,x= x,color="imcem"),size=0.5) +
  scale_color_manual(values = c(
    'saem' = 'darkblue','mcem' = 'red','imcem' = 'black')) +
  labs(color = 'Algo')
