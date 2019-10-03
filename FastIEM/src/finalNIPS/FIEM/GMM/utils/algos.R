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
    # if (k %% n==0)
    # {
    #   print('EM')
    #   print(k)
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


mixt.iem.seq <- function(x, theta0, K,nbr)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  tau <- compute.tau(x,theta0)
  theta<-theta0
  # tau.old <- compute.tau(x[1],theta0)
  s <- compute.stat(x,tau)
  l <- rep(1:n,K/n)
  i <- 1
  for (k in 1:K)
  {

    # if (k %% n==0)
    # {
    #   print('IEM')
    #   print(k)
    # }
    # i <- sample(1:n, 1)
    #Update the conditional expectation for the chosen datum
    oldtau <- tau[l[i],]
    tau[l[i],] <- compute.tau(x[l[i]],theta)
    
    #Update the statistics 
    s$s1 <- s$s1 + tau[l[i],] - oldtau
    s$s2 <- s$s2 + x[l[i]]*(tau[l[i],] - oldtau)
    # s <- compute.stat(x,tau)
    
    #M-step
    theta$mu<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)

    i <- i +nbr
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
  s <- compute.stat(x,tau)

  for (k in 1:K)
  {

    # if (k %% n==0)
    # {
    #   print('IEM')
    #   print(k)
    # }
    i <- sample(1:n, 1)
    #Update the conditional expectation for the chosen datum
    oldtau <- tau[i,]
    tau[i,] <- compute.tau(x[i],theta)
    
    #Update the statistics 
    s$s1 <- s$s1 + tau[i,] - oldtau
    s$s2 <- s$s2 + x[i]*(tau[i,] - oldtau)
    # s <- compute.stat(x,tau)
    
    #M-step
    theta$mu<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)
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
  
  for (k in 1:K)
  {

    # if (k %% n==0)
    # {
    #   print('OEM')
    #   print(k)
    # }
    i <- sample(1:n, 1)
    tau.indiv.new <- compute.tau(x[i],theta)
    s.indiv.new <- x[i]*tau.indiv.new

    #Update statistic
    s$s1 <- s$s1 + rho[k]*(tau.indiv.new  - s$s1)
    s$s2 <- s$s2 + rho[k]*(s.indiv.new  - s$s2)
    
    #M-step
    theta$mu <- step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)

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
  

  theta.e.0 <- theta
  tau.e.0 <- compute.tau(x,theta.e.0)
  s.e.0 <- x%*%tau.e.0


  for (k in 1:K)
  {

    if (k%%(n/nbr) == 0)
    { 
      # print('OEMVR')
      # print(k)
      theta.e.0 <- theta
      tau.e.0 <- compute.tau(x,theta.e.0)
      s.e.0 <- x%*%tau.e.0
    }
    i <- sample(1:n, 1)
    tau.indiv.new <- compute.tau(x[i],theta)
    tau.indiv.e.0 <- compute.tau(x[i],theta.e.0)

    #Update statistics
    s$s1 <- (1-rho)*s$s1 + rho*((tau.indiv.new - tau.indiv.e.0)*n + colSums(tau.e.0))
    s$s2 <- (1-rho)*s$s2 + rho*((x[i]*tau.indiv.new - x[i]*tau.indiv.e.0)*n + s.e.0)

    #M-step
    theta$mu <- step.M(s,n)
    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)

  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.saga <- function(x, theta0, K,nbr, rho.saga)
{
   G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  
  #Init
  tau <- compute.tau(x,theta)
  s <- v <- h <- compute.stat(x,tau)
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
    # if (k %% n==0)
    # {
    #   print('SAGA')
    #   print(k)
    # }

    newtau.i<- compute.tau(x[li[i]],theta)
    oldtau.i<- compute.tau(x[li[i]],alphas[[li[i]]])


    v$s1 <- h$s1 + (newtau.i - oldtau.i)*n
    v$s2 <- h$s2 + (x[li[i]]*newtau.i - x[li[i]]*oldtau.i)*n
    
  
    s$s1 <- (1-rho.saga)*s$s1 + rho.saga*v$s1
    s$s2 <- (1-rho.saga)*s$s2 + rho.saga*v$s2

    oldtheta <- theta
    theta$mu<-step.M(s,n)

    theta.est[k+1,] <- c(k, theta0$p, theta$mu, theta0$sigma)

    oldalpha.j <- alphas[[lj[j]]]
    alphas[[lj[j]]] <- oldtheta
    newtau.j<- compute.tau(x[lj[j]],oldtheta)
    oldtau.j<- compute.tau(x[lj[j]],oldalpha.j)
    # tau[lj[j],] <- newtau.i - oldtau.i
    h$s1 <- h$s1 + (newtau.j - oldtau.j)
    h$s2 <- h$s2 + (x[lj[j]]*newtau.j - x[lj[j]]*oldtau.j)

    i <- i+nbr
    j <- j+nbr
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}





