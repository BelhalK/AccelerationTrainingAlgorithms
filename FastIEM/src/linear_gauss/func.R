require(ggplot2)
require(gridExtra)
require(reshape2)


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

mixt.ident3 <- function(df)
{
  G <- (ncol(df)-1)/3
  K <- nrow(df)
  mu.final <- as.numeric(as.character(df[K,(G+2):(2*G+1)]))
  ind <- sort.int(mu.final, index.return=TRUE)$ix
  df[,2:(G+1)] <- df[,(G-2+ind)]
  df[,(G+2):(2*G+1)] <- df[,(2*G-2+ind)]
  df[,(2*G+2):(3*G+1)] <- df[,(3*G-2+ind)]
  return(df)
}


graphConvMC_new <- function(df, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$rep <- as.factor(df$rep)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +
      xlab("iteration") + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}

graphConvMC2_new <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df,aes(colour=df$algo ))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend) +
      xlab("iteration") + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
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


compute.tau<-function(x,theta, alph)
{
  n<-length(x)
  tau<-matrix(NA,n,1)
  for (i in (1:n))
  {
    tau[i,]<-alph*theta$mu[1]+(1-alph)*x[i]
  }
  return(tau)
}
  

compute.stat<-function(x,Z)
{
  M<-dim(Z)[3]
  if (is.na(M))  
  {
    M <- 1
    dim(Z) <- c(dim(Z),1)
  }
  s1 <- 0
  for (m in 1:M)
  {
    Z.m <- Z[,,m]
    s1 <- s1 + sum(Z.m)
  }
  s <-list(s1=s1/M)
  return(s)
}


step.E<-function(x,theta, alph)
{
  tau <- compute.tau(x,theta, alph)
  s <- compute.stat(x,tau)
  return(s)
}

step.M<-function(s,n)
{
  mu<-s$s1/n
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