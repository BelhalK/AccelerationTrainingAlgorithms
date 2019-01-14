require(ggplot2)
require(gridExtra)
require(reshape2)

mixt.simulate <-function(n,mu,sigma)
{
  x<-NULL
  x<-rnorm(n,mu[1],sigma[1]+sigma[2])
  return(x)
}



#-------------------------------------
mixt.em <- function(x, theta0, K, alph)
{
  G<-1
  col.names <- c("iteration", paste0("mu",1:G))
  theta.est <- matrix(NA,K+1,2)
  theta.est[1,] <- c(0, theta0$mu)
  theta<-theta0
  for (k in 1:K)
  {
    s<-step.E(x,theta,alph)
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k,theta$mu)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.iem <- function(x, theta0, K, alph,nbr)
{
  G<-1
  col.names <- c("iteration", paste0("mu",1:G))
  theta.est <- matrix(NA,K+1,2)
  theta.est[1,] <- c(0, theta0$mu)
  theta<-theta0
  tau<-compute.tau(x,theta,alph)
  n<-length(x)
  l <- NULL
  # for (j in 1:(K/n))
  # {
  #   l <- c(l, sample(1:n,n))
  # }
  # l <- sample(1:n,K,replace = TRUE)
  # l <- rep(1:n,K/n)
  l <- rep(sample(1:n,n), K/n)
  i <- 1:nbr
  for (k in 1:K)
  {
    if (k%%(n/nbr) == 1)
    { 
      # l<-sample(1:n,n)
      # l<-1:n
      i<-1:nbr
    }

    # i<-sample(1:n,nbr)
    tau[l[i],] <- compute.tau(x[l[i]],theta,alph)
    # tau[i,] <- compute.tau(x[i],theta,alph)
    s <- compute.stat(x,tau)
    i <- i+nbr
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k,theta$mu)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

mixt.oem <- function(x, theta0, K, alph,nbr)
{
  G<-1
  kiter = 1:K
  rho = 3/(kiter+10)
  col.names <- c("iteration", paste0("mu",1:G))
  theta.est <- matrix(NA,K+1,2)
  theta.est[1,] <- c(0, theta0$mu)
  theta<-theta0
  tau<-compute.tau(x,theta,alph)
  s <- compute.stat(x,tau)
  n<-length(x)
  l <- NULL
  # for (j in 1:(K/n))
  # {
  #   l <- c(l, sample(1:n,n))
  # }
  # l <- sample(1:n,K,replace = TRUE)
  # l <- rep(1:n,K/n)
  l <- rep(sample(1:n,n), K/n)
  i <- 1:nbr
  for (k in 1:K)
  {
    if (k%%(n/nbr) == 1)
    { 
      # l<-sample(1:n,n)
      # l<-1:n
      i<-1:nbr
    }

    # i<-sample(1:n,nbr)
    s.old <- compute.stat(x,tau)
    tau[l[i],] <- compute.tau(x[l[i]],theta,alph)
    # s <- compute.stat(x,tau)
    s$s1 <- s.old$s1 + rho[k]*(tau[l[i],] - s.old$s1)
    i <- i+nbr
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k,theta$mu)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.oemvr <- function(x, theta0, K, alph,nbr)
{
  G<-1
  rho =0.1
  col.names <- c("iteration", paste0("mu",1:G))
  theta.est <- matrix(NA,K+1,2)
  theta.est[1,] <- c(0, theta0$mu)
  theta<-theta0
  tau<-compute.tau(x,theta,alph)
  tau.old <- tau[1,]
  s <- compute.stat(x,tau)
  s.old.init <- s
  n<-length(x)
  
  l <- NULL
  l <- rep(sample(1:n,n), K/n)
  i <- 1:nbr
  
  for (k in 1:K)
  {
    if (k%%(n/nbr) == 1)
    { 
      tau.old <- compute.tau(x[l[i]],theta,alph)
      s.old.init <- s
      i<-1:nbr
    }

    # i<-sample(1:n,nbr)
    s.old <- compute.stat(x,tau)
    tau[l[i],] <- compute.tau(x[l[i]],theta,alph)
    # s <- compute.stat(x,tau)
    s$s1 <- s.old$s1 + rho*(tau[l[i],] - tau.old + s.old.init$s1 - s.old$s1)
    i <- i+nbr
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k,theta$mu)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}
