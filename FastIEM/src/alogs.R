require(ggplot2)
require(gridExtra)
require(reshape2)

mixt.em <- function(x, theta0, K,A,B,sigma,omega,id)
{
  G<-length(theta0$mu)
  n <- length(unique(id))
  col.names <- c("iteration", paste0("mu",1:G))
  theta.est <- matrix(NA,K+1,3)
  theta.est[1,] <- c(0, theta0$mu)
  theta<-theta0
  for (k in 1:K)
  {
    s<-step.E(x,theta,A,B,sigma,omega,id)
    theta<-step.M(s,n,A,B,sigma,x,id)
    theta.est[k+1,] <- c(k,theta$mu)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.iem <- function(x, theta0, K,A,B,sigma,omega,id,nbr)
{
  G<-length(theta0$mu)
  n <- length(unique(id))
  ni <- nrow(A)
  col.names <- c("iteration", paste0("mu",1:G))
  theta.est <- matrix(NA,K+1,3)
  theta.est[1,] <- c(0, theta0$mu)
  theta<-theta0

  tau <- matrix(0,nrow=n,ncol=nrow(omega))
  somega <- solve(omega)
  ssigma <- solve(sigma)
  gamma <- solve(t(B)%*%ssigma%*%B+somega)
  const <- gamma%*%t(B)%*%ssigma
  for (j in (1:n))
  {  
    tau[j,]<-const%*%(x[j,] - A%*%theta$mu)
  }

  l <- NULL
  l<-rep(sample(1:n,n),K)
  i <- 1:nbr

  for (k in 1:K)
  {
    if (nbr==n){
      l<-1:n
      i<-1:n
    } else{
      if (k%%(n/nbr) == 0)
      { 
        # l<-sample(1:n,n)
        # l<-1:n
        i<-1:nbr
      }
    }
    for (j in l[i])
    {  
      tau[j,]<-const%*%(x[j,] - A%*%theta$mu)
    }
    s<-colMeans(tau)
    theta<-step.M(s,n,A,B,sigma,x,id)
    theta.est[k+1,] <- c(k,theta$mu)
    i <- i+nbr
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}





# mixt.iem_seq <- function(x, theta0, K,A,B,sigma,omega,id,nbr)
# {
#   G<-length(theta0$mu)
#   n <- length(unique(id))
#   col.names <- c("iteration", paste0("mu",1:G))
#   theta.est <- matrix(NA,K+1,3)
#   theta.est[1,] <- c(0, theta0$mu)
#   theta<-theta0
#   tau <- compute.tau(x,theta,A,B,sigma,omega,id)
#   l <- NULL
#   # for (j in 1:(K/n))
#   # {
#   #   l <- c(l, sample(1:n,n))
#   # }
#   # l <- sample(1:n,K,replace = TRUE)
#   # l <- rep(1:n,K/n)
#   l <- rep(sample(1:n,n), K/n)
#   i <- 1:nbr
#   for (k in 1:K)
#   {
#     if (k%%(n/nbr) == 1)
#     { 
#       # l<-sample(1:n,n)
#       # l<-1:n
#       i<-1:nbr
#     }
#     # i<-sample(1:n,nbr)
#     tau[l[i],] <- compute.tau_incremental(x[which(id==l[i])],theta,A,B,sigma,omega,l[i])
#     # tau[i,] <- compute.tau(x[i],theta,alph)
#     s <- compute.stat(x,tau)
#     i <- i+nbr
#     theta<-step.M(s,n,A,B,sigma,x,id)
#     theta.est[k+1,] <- c(k,theta$mu)
#   }
  
#   df <- as.data.frame(theta.est)
#   names(df) <- col.names
#   return(df)
# }

# mixt.iem_periter <- function(x, theta0, K,A,B,sigma,omega,id,nbr)
# {
#   G<-length(theta0$mu)
#   n <- length(unique(id))
#   col.names <- c("iteration", paste0("mu",1:G))
#   theta.est <- matrix(NA,K+1,3)
#   theta.est[1,] <- c(0, theta0$mu)
#   theta<-theta0
#   tau <- compute.tau(x,theta,A,B,sigma,omega,id)
#   l <- NULL
#   for (k in 1:K)
#   {

#     i<-sample(1:n,nbr)
#     tau[i,] <- compute.tau_incremental(x[which(id==i)],theta,A,B,sigma,omega,l[i])
#     s <- compute.stat(x,tau)
#     i <- i+nbr
#     theta<-step.M(s,n,A,B,sigma,x,id)
#     theta.est[k+1,] <- c(k,theta$mu)
#   }
  
#   df <- as.data.frame(theta.est)
#   names(df) <- col.names
#   return(df)
# }


mixt.iem_perpass <- function(x, theta0, K, alph,nbr)
{
  G<-1
  col.names <- c("iteration", paste0("mu",1:G))
  theta.est <- matrix(NA,K+1,2)
  theta.est[1,] <- c(0, theta0$mu)
  theta<-theta0
  tau<-compute.tau(x,theta,alph)
  n<-length(x)
  l <- NULL
  for (k in 1:K)
  {
    if (k%%(n/nbr) == 1)
    { 
      l<-sample(1:n,n)
      # l<-1:n
      i<-1:nbr
    }

    tau[l[i],] <- compute.tau(x[l[i]],theta,alph)
    s <- compute.stat(x,tau)
    i <- i+nbr
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k,theta$mu)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

#-------------------------------------
mixt.piem <- function(x, theta0, K,alph)
{
  G<-1
  col.names <- c("iteration", paste0("mu",1:G))
  theta.est <- matrix(NA,K+1,2)
  theta.est[1,] <- c(0, theta0$mu)
  theta<-theta0
  tau<-compute.tau(x,theta,alph)
  liste <- matrix(NA,K,1)
  tau_old <- tau
  for (k in 1:K)
  {
    diff <- matrix(NA,n,1)
    tau_test <- compute.tau(x,theta,alph)
    diff <- (tau_test[,1] - tau_old[,1])^2
    i <- which.max(diff)
    liste[k] <- i
    tau[i,] <- compute.tau(x[i],theta,alph)
    tau_old <- tau
    s <- compute.stat(x,tau)
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$mu)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}
#-------------------------------------
