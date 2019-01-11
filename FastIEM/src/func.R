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

graphConvMC <- function(df, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$rep <- as.factor(df$rep)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +
      xlab("iteration") + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
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
    s3 <- s3 + (x^2) %*% Z.m 
  }
  s <-list(s1=s1/M,s2=as.vector(s2/M),s3=as.vector(s3/M))
  return(s)
}

compute.stat_iem<-function(x,Z, tau.new, tau.old, i)
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
    s1 <- s1 + colSums(Z.m) + (tau.new- tau.old)
    s2 <- s2 + x %*% Z.m + x[i]*(tau.new- tau.old)
    s3 <- s3 + (x^2) %*% Z.m + x[i]^2*(tau.new- tau.old)
  }
  s <-list(s1=s1/M,s2=as.vector(s2/M),s3=as.vector(s3/M))
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
  p<-s$s1/n
  mu<-s$s2/s$s1
  sigma<-sqrt(s$s3/s$s1-(mu)^2)
  theta<-list(p=p,mu=mu,sigma=sigma)
  return(theta)
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

