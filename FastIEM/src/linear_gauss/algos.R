require(ggplot2)
require(gridExtra)
require(reshape2)

mixt.simulate <-function(n,weight,mu,sigma)
{
  G <- length(mu)
  Z <- sample(1:G, n, prob=weight, replace=T)
  x<-NULL
  for (g in 1:G)
  {
    x<-c(x,rnorm(length(which(Z==g)),mu[g],sigma[g]))
  }
  return(x)
}



#-------------------------------------

mixt.em <- function(x, theta0, K)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  
  theta<-theta0
  for (k in 1:K)
  {
    # if (k %% 1000==0)
    # {
    #   print(k)
    # }
    s<-step.E(x,theta)
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

mixt.iem <- function(x, theta0, K,nbr)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  tau <- compute.tau(x,theta0)
  theta<-theta0
  tau.old <- compute.tau(x[1],theta0)
  s <- compute.stat_iem(x,tau, tau.old, tau.old,1)

  l <- rep(sample(1:n,n), K/n)
  i <- 1:nbr
  for (k in 1:K)
  {
    print(k)
    if (k%%(n/nbr) == 1)
    { 
      # l<-sample(1:n,n)
      # l<-1:n
      i<-1:nbr
    }
    tau.new <- compute.tau(x[i],theta)
    s <- compute.stat_iem(x,tau, tau.new,tau.old, i)
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
    i = i+nbr
    tau.old <- tau.new
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.oem <- function(x, theta0, K, alph,nbr)
{
  G<-1
  kiter = 1:K
  rho = 10/(kiter+10)
  col.names <- c("iteration", paste0("mu",1:G))
  theta.est <- matrix(NA,K+1,2)
  theta.est[1,] <- c(0, theta0$mu)
  theta<-theta0
  tau<-compute.tau(x,theta,alph)
  s<-step.E(x,theta,alph)
  theta<-step.M(s,n)
  n<-length(x)
  l <- NULL
  l <- rep(sample(1:n,n), K/n)
  i <- 1:nbr
  for (k in 1:K)
  {
    if (k%%(n/nbr) == 1)
    { 
      i<-1:nbr
    }
    tau.i <- compute.tau(x[l[i]],theta,alph)
    s$s1 <- s$s1 + rho[k]*(tau.i - s$s1)
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
  kiter = 1:K
  rho =0.01
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
    s.old <- s
    if (k%%(n/nbr) == 1)
    { 
      i<-1:nbr
      tau.old <- compute.tau(x[l[i]],theta,alph)
      s.old.init <- s.old
    }

    tau[l[i],] <- compute.tau(x[l[i]],theta,alph)
    s$s1 <- s.old$s1 + rho*(tau[l[i],] - tau.old + s.old.init$s1 - s.old$s1)
    i <- i+nbr
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k,theta$mu)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}
