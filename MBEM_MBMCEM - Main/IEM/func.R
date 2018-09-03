require(ggplot2)
require(gridExtra)
require(reshape2)


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


step.E_incremental<-function(x,theta, A,B,sigma,omega,id)
{
   n<-length(unique(id))
  tau <- matrix(0,nrow=n,ncol=nrow(omega))
  somega <- solve(omega)
  ssigma <- solve(sigma)
  gamma <- solve(t(B)%*%ssigma%*%B+somega)
  const <- gamma%*%t(B)%*%ssigma
  for (i in (1:n))
  {  
    tau[i,]<-const%*%(x[i,] - A%*%theta$mu)
  }
  return(colMeans(tau))
}


step.E<-function(x,theta, A,B,sigma,omega,id)
{
   n<-length(unique(id))
  tau <- matrix(0,nrow=n,ncol=nrow(omega))
  somega <- solve(omega)
  ssigma <- solve(sigma)
  gamma <- solve(t(B)%*%ssigma%*%B+somega)
  const <- gamma%*%t(B)%*%ssigma
  for (i in (1:n))
  {  
    tau[i,]<-const%*%(x[i,] - A%*%theta$mu)
  }
  return(colMeans(tau))
}

step.M<-function(s,n,A,B,sigma,x,id)
{ 
  obs <- colMeans(x)
  ssigma <- solve(sigma)
  int <- solve(t(A)%*%ssigma%*%A)%*%t(A)%*%ssigma
  mu<-int%*%(obs - B%*%s)
  theta<-list(mu=mu)
  return(theta)
}

#simulation

step.S<-function(x,theta,M,alph,gamm)
{
  n<-length(x)
  G<-1
  Z<-array(NA,c(n,G,M))
  tau<-compute.tau(x,theta,alph)
  for (i in 1:n)
    Z[i,,]<-rnorm(M,tau[i,],gamm)
  return(Z)
}


#stepStochasticApproximation

step.SA <-function(x,Z,s.old,gamma)
{
  S<-compute.stat(x,Z)
  s11<-s.old$s1+gamma*(S$s1-s.old$s1)
  s.new<-list(s1=s11)
  return(s.new)
}

step.SAll <-function(x,Z,s.old,gamma)
{
  # S<-compute.stat(x,Z)
  s11<-s.old+gamma*(Z-s.old)
  # s.new<-list(s1=s11)
  return(s11)
}

step.SAmeanfield <-function(x,tau,s.old,gamma)
{
  n<-length(x)
  s11 <- matrix(NA,n,1)
  for (i in 1:n)
  {
    s11[i,]<-1/n*tau[i,]-s.old[i,]
  }
  # s.new<-list(s1=s11)
  return(s11)
}


step.Mh<-function(h,n)
{
  s1 <- 0
  s1 <- sum(h)
  mu<-s1/n
  theta<-list(mu=mu)
  return(theta)
}

