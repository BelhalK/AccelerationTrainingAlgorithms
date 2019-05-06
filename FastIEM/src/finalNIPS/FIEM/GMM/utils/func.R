

mixt.ident <- function(df)
{
  G <- (ncol(df)-1)/3
  K <- nrow(df)
  mu.final <- as.numeric(as.character(df[K,(G+2):(2*G+1)]))
  ind <- sort.int(mu.final, index.return=TRUE)$ix
  df[,2:(G+1)] <- df[,(G-1+ind)]
  df[,(G+2):(2*G+1)] <- df[,(2*G-1+ind)]
  df[,(2*G+2):(3*G+1)] <- df[,(3*G-1+ind)]
  return(df)
}


logLikelihood <- function(x,df,sigma)
{
  K <- dim(df)[1]
  G <- 1
  ll <- NULL
  for (k in (1:K))
  {
    lk <- 0
    mu.k <- df[k,2]
    lk <- dnorm(x, mean=mu.k, sd=sigma[1]+sigma[2],log=TRUE)
    ll <- c(ll, sum(lk))
  }
  df$deviance <- -2*ll
  return(df)
}



compute.tau<-function(x,theta)
{
  n<-length(x)
  G<-length(theta$p)
  tau<-matrix(NA,n,G)
  for (g in 1:G)
    for (i in (1:n))
    {
      tau[i,g]<-theta$p[g]*dnorm(x[i],theta$mu[g],theta$sigma[g])
    }
  # tau<-prop.table(tau,1)
  tau=tau/matrix(rep(rowSums(tau),G),nrow=n)
  return(tau)
}



compute.stat<-function(x,Z)
{
  G<-dim(Z)[2]
  M<-dim(Z)[3]
  if (is.na(M))  
  {
    M <- 1
    dim(Z) <- c(dim(Z),1)
  }
  s1 <- 0
  s2 <- 0
  s3 <- 0
  for (m in 1:M)
  {
    Z.m <- Z[,,m]
    s1 <- s1 + colSums(Z.m) 
    s2 <- s2 + x %*% Z.m 
    # s3 <- s3 + (x^2) %*% Z.m 
  }
  s <-list(s1=s1/M,s2=as.vector(s2/M))
  return(s)
}

step.E<-function(x,theta)
{
  tau <- compute.tau(x,theta)
  s <- compute.stat(x,tau)
  return(s)
}

step.M<-function(s,n)
{
  # p<-s$s1/n
  mu<-s$s2/s$s1
  # sigma<-sqrt(s$s3/s$s1-(mu)^2)
  # theta<-list(mu=mu)
  return(mu)
}



step.S<-function(x,theta,M)
{
  n<-length(x)
  G<-length(theta$p)
  Z<-array(NA,c(n,G,M))
  tau<-compute.tau(x,theta)
  test <- F
  while (test==F)
  { 
    for (i in 1:n)
      Z[i,,]<-rmultinom(n=M,size=1,prob=tau[i,])
    test <- (min(colSums(Z[,,1]))>1) 
  }
  return(Z)
}

