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
    theta$mu<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)
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
  # tau.old <- compute.tau(x[1],theta0)
  # s <- compute.stat_iem(x,tau, tau.old,1)

  l <- rep(sample(1:n,n), K/n)
  i <- 1:nbr
  for (k in 1:K)
  {

    if (k%%(n/nbr) == 1)
    { 
      i<-1:nbr
    }
    
    tau[l[i],] <- compute.tau(x[l[i]],theta)
    s <- compute.stat(x,tau)
    theta$mu<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)
    i = i+nbr
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.oem <- function(x, theta0, K,nbr)
{
  G<-length(mu)
  kiter = 1:K
  rho = 3/(kiter+10)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  tau<-compute.tau(x,theta)
  s<-compute.stat(x,tau)
  # s<-step.E(x,theta)
  # theta$mu<-step.M(s,n)
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

    tau.new <- compute.tau(x[l[i]],theta)
    s <- compute.stat_oem(x,tau, tau.new, l[i],rho[k])

    i <- i+nbr
    theta$mu<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)
    # tau[l[i],] <- tau.new
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

mixt.oemvr <- function(x, theta0, K,nbr,rho)
{
   G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  
  tau <- compute.tau(x,theta)
  tau.old <- tau[1,]
  
  s<-compute.stat(x,tau)
  s.old.init <- s
  
  l <- NULL
  l <- rep(sample(1:n,n), K/n)
  i <- 1:nbr
  
  for (k in 1:K)
  {
    if (k%%(n/nbr) == 1)
    { 
      i<-1:nbr
      tau.old <- compute.tau(x[l[i]],theta)
      s.old.init <- compute.stat(x,tau)
    }
    
    tau.new <- compute.tau(x[l[i]],theta)
    s <- compute.stat_oemvr(x,tau, tau.new,s.old.init,tau.old, l[i],rho)
    tau[l[i],] <-tau.new

    # s.old<-compute.stat(x,tau)
    # tau[l[i],] <- compute.tau(x[l[i]],theta)
    # s1 <- colSums(tau)
    # s$s2 <- s.old$s2 + rho*(x[i]%*%tau[l[i],] - x[l[i]]%*%tau.old + s.old.init$s2 - s.old$s2)
    # s$s2 <- s$s2/s1

    i <- i+nbr
    theta$mu<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


