library(MASS)
library(mvtnorm)
source("MCEM_mulvar_fcts.R")

N <- 1000  # number of subjects
p <- 5     # number of explanatory variables

mu.star <- 1:p  #rep(0,p)  # mean of the explanatory variables
sd <- 1:p # rep(1,p) # standard deviations

C <- matrix(c(   # correlation matrix
  1,   0.8, 0,   0,   0,
  0.8, 1,   0,   0,   0,
  0,   0,   1,   0.3, 0.6,
  0,   0,   0.3, 1,   0.7,
  0,   0,   0.6, 0.7, 1
), nrow=p)

Sigma.star <- diag(sd)%*%C%*%diag(sd) # variance-covariance matrix of the explanatory variables

beta.star <- c(0.5, -0.3, 1, 0, -0.6) # coefficients
beta0.star <- -0.2  # intercept

p.miss <- 0.25  # fraction of mising data

nbsim = 2
#e.obs = e.predx = e.predxy = e.saem = f.predx= f.predxy = f.saem = NULL
EST.comp = EST.cc = EST.mcem.ar = matrix(0, nbsim,length(beta.star)+1)
BIAS.comp = BIAS.cc = BIAS.mcem.ar = rep(0, nbsim)

#1 simu
X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) + matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
y <- as.numeric(runif(N)<p1)

data.complete <- data.frame(y=y,X.complete)
model.complete <- glm(y ~.,family=binomial(link='logit'),data=data.complete)
beta0.complete <- model.complete$coefficients[1]
beta.complete <- model.complete$coefficients[2:(p+1)]
bias.complete  = sum((c(beta0.complete ,beta.complete)-c(beta0.star,beta.star))^2)

# ------- generating missing data
X.obs <- X.complete
X.obs[runif(N*p)<p.miss] <- NA



#simus
for (NB in 1:nbsim){
  cat(sprintf('simulation = %i ', NB,'/n'))

  # ----- complete data simulation
  X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) + matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
  p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
  y <- as.numeric(runif(N)<p1)

  data.complete <- data.frame(y=y,X.complete)
  model.complete <- glm(y ~.,family=binomial(link='logit'),data=data.complete)
  beta0.complete <- model.complete$coefficients[1]
  beta.complete <- model.complete$coefficients[2:(p+1)]
  bias.complete  = sum((c(beta0.complete ,beta.complete)-c(beta0.star,beta.star))^2)

  # ------- generating missing data
  X.obs <- X.complete
  X.obs[runif(N*p)<p.miss] <- NA

  # ------- estimation ignoring the missing data
  data.obs <- data.frame(y=y,X.obs)
  model.obs <- glm(y ~.,family=binomial(link='logit'),data=data.obs)
  beta0.cc <- model.obs$coefficients[1]
  beta.cc <- model.obs$coefficients[2:(p+1)]
  bias.cc = sum((c(beta0.cc,beta.cc)-c(beta0.star,beta.star))^2)


  #list.est = mcem_ar(X= X.obs,Y= y , maxruns=1000,tol_em=1e-5)
  #ptm <- proc.time()
  list.est = mcem_ar(X= X.obs,Y= y, maxruns=1000, tol_em=1e-3)
  #time_mcem_ar=proc.time() - ptm
  bias.est = sum((list.est$beta.est-c(beta0.star,beta.star))^2)

  EST.comp[NB,] = c(beta0.complete,beta.complete)
  EST.cc[NB,] = c(beta0.cc,beta.cc)
  EST.mcem.ar[NB,] = list.est$beta.est
  BIAS.comp[NB] = bias.complete
  BIAS.cc[NB] = bias.cc
  BIAS.mcem.ar[NB] = bias.est
}
