require(ggplot2)
require(gridExtra)
require(reshape2)
library(rlist)
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
      # print(k)
    # }
    
    #Update the statistics
    s<-step.E(x,theta)

    #M-step
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
      # print(k)
    }
    
    #Update the conditional expectation for the chosen datum
    tau[l[i],] <- compute.tau(x[l[i]],theta)
    
    #Update the statistics
    s <- compute.stat(x,tau)
    
    #M-step
    theta$mu<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)
    i = i+nbr
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.oem <- function(x, theta0, K,nbr,rho)
{
   G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0

  #Init
  tau <- compute.tau(x,theta)
  s<-compute.stat(x,tau)
  
  l <- NULL
  l <- rep(sample(1:n,n), K/n)
  i <- 1:nbr
  
  for (k in 1:K)
  {
    if (k%%(n/nbr) == 1)
    { 
      i<-1:nbr
      # print(k)
    }
    tau.indiv.new <- compute.tau(x[l[i]],theta)
    s.indiv.new <- x[l[i]]*tau.indiv.new

    #Update statistic
    s$s1 <- s$s1 + rho[k]*(tau.indiv.new  - s$s1)
    s$s2 <- s$s2 + rho[k]*(s.indiv.new  - s$s2)
    
    #M-step
    theta$mu <- step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)

    #Update index
    i <- i+nbr

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
  
  #Init
  tau <- compute.tau(x,theta)
  s<-compute.stat(x,tau)
  
  l <- NULL
  l <- rep(sample(1:n,n), K/n)
  i <- 1:nbr
  
  for (k in 1:K)
  {
    if (k%%(n/nbr) == 1)
    { 
      i<-1:nbr
      theta.e.0 <- theta
      # print(k)
    }
    tau.indiv.new <- compute.tau(x[l[i]],theta)
    s.indiv.new <- x[l[i]]*tau.indiv.new

    tau.indiv.e.0 <- compute.tau(x[l[i]],theta.e.0)
    s.indiv.e.0 <- x[l[i]]*tau.indiv.e.0

    tau.e.0 <- compute.tau(x,theta.e.0)
    s.e.0 <- x%*%tau.e.0

    #Update statistic
    s$s1 <- s$s1 + rho*(tau.indiv.new - tau.indiv.e.0 + colSums(tau.e.0) - s$s1)
    s$s2 <- s$s2 + rho*(s.indiv.new - s.indiv.e.0 + s.e.0 - s$s2)

    #M-step
    theta$mu <- step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)

    #Update index
    i <- i+nbr
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.saga <- function(x, theta0, K,nbr)
{
   G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  
  #Init
  tau <- compute.tau(x,theta)
  s <- compute.stat(x,tau)
  v <- compute.stat(x,tau)
  n<-length(x)
  li <- NULL
  alphas <- rep(list(theta0),n)
  # l <- sample(1:n,K,replace = TRUE)
  # l <- rep(1:n,K/n)
  li <- rep(sample(1:n,n), K/n)
  lj <- NULL
  for (index in 1:(K/n)){
    lj <- list.append(lj, sample(li[(1+(index-1)*n):(index*n)]))
  }
  i <- 1:nbr
  j <- 1:nbr
  
  for (k in 1:K)
  {
    if (k%%(n/nbr) == 1)
    { 
      i<-1:nbr
      j<-1:nbr
    }
    newtau.i<- compute.tau(x[li[i]],theta)
    oldtau.i<- compute.tau(x[li[i]],alphas[[li[i]]])
    # tau[li[i],] <- (newtau.i - oldtau.i)*n

    v$s1 <- s$s1 + (newtau.i - oldtau.i)*n
    v$s2 <- s$s2 + (x[li[i]]*newtau.i - x[li[i]]*oldtau.i)*n
    
    oldtheta <- theta
    theta$mu<-step.M(v,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)

    oldalpha.j <- alphas[[lj[j]]]
    alphas[[lj[j]]] <- oldtheta
    newtau.j<- compute.tau(x[lj[j]],oldtheta)
    oldtau.j<- compute.tau(x[lj[j]],oldalpha.j)
    # tau[lj[j],] <- newtau.i - oldtau.i
    s$s1 <- s$s1 + (newtau.j - oldtau.j)
    s$s2 <- s$s2 + (x[lj[j]]*newtau.j - x[lj[j]]*oldtau.j)

    i <- i+nbr
    j <- j+nbr
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}





