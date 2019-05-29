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
  df$sim <- as.factor(df$sim)
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
  do.call("grid.arrange", c(graf, ncol=3, top=title))
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
      xlab("iteration") +scale_x_log10()+ ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
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

graphConvMC2 <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df,aes(colour=df$algo ))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend) +
      xlab("iteration")+scale_x_log10() + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}

graphConvMC3 <- function(df, title=NULL, ylim=NULL,legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$rep <- as.factor(df$rep)
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-3)
  o <- c(0, 1, 2, 3, 4, 5, 6)
  for (j in (2:(ncol(df)-2)))
  {
    grafj <- ggplot(df,aes(colour=df$rep, linetype=df$algo))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend) + 
    xlab("iteration") + ylab(names(df[j]))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}


graphSummary <- function(D)
{
  n <- length(D)
  G <- (ncol(df)-2)/3
  df$rep <- as.factor(df$rep)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +
      xlab("iteration") + ylab(names(df[j])) + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}

graphConvergence <- function(df, title=NULL, ylim=NULL)
{
  G <- (dim(df)[2]-1)/3
  df1<-melt(df[,c(1,2,3)],id.vars="iteration",value.name="proportion",variable.name="group")
  graf1 <- ggplot(df1)+geom_line(aes(iteration,proportion,color=group)) + theme(legend.position="none")  
  if (!is.null(ylim))  graf1 <- graf1 + ylim(ylim[1]*c(-1,1))
  df2<-melt(df[,c(1,4,5)],id.vars="iteration",value.name="mean",variable.name="group")
  graf2 <- ggplot(df2)+geom_line(aes(iteration,mean,color=group)) + theme(legend.position="none") 
  if (!is.null(ylim))  graf2 <- graf2 + ylim(ylim[2]*c(-1,1))
  df3<-melt(df[,c(1,6,7)],id.vars="iteration",value.name="standard.deviation",variable.name="group")
  graf3 <- ggplot(df3) + geom_line(aes(iteration,standard.deviation,color=group)) + theme(legend.position="none") 
  if (!is.null(ylim))  graf3 <- graf3 + ylim(ylim[3]*c(-1,1))
  graf4 <- ggplot(df)+geom_line(aes(iteration,deviance))  
  grid.arrange(graf1,graf2,graf3,graf4,nrow=2, top=title)
}

logLikelihood <- function(x,df)
{
  K <- dim(df)[1]
  G <- (dim(df)[2]-1)/3
  ll <- NULL
  for (k in (1:K))
  {
    lk <- 0
    for (g in (1:G))
    {
      pi.k <- df[k,g+1]
      mu.k <- df[k,G+g+1]
      sigma.k <- df[k,2*G+g+1]
      lk <- lk + pi.k * dnorm(x, mean=mu.k, sd=sigma.k)
    }
    ll <- c(ll, sum(log(lk)))
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
    tau[,g]<-theta$p[g]*dnorm(x,theta$mu[g],theta$sigma[g])
  # tau<-prop.table(tau,1)
  tau=tau/matrix(rep(rowSums(tau),G),nrow=n)
  return(tau)
}

compute.tau_replace<-function(x,theta,nb_r)
{
  n<-length(x)
  G<-length(theta$p)
  tau<-matrix(NA,nb_r,G)
  for (g in 1:G)
    tau[,g]<-theta$p[g]*dnorm(x,theta$mu[g],theta$sigma[g])
  # tau<-prop.table(tau,1)
  tau=tau/matrix(rep(rowSums(tau),G),nrow=nb_r)
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

#simulation
#stepSimulation
step.S_replace<-function(x,theta,M)
{
  n<-length(x)
  G<-length(theta$p)
  Z<-array(NA,c(n,G,M))
  tau<-compute.tau(x,theta)
  for (i in 1:n)
    Z[i,,]<-rmultinom(n=M,size=1,prob=tau[i,])
  return(Z)
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

step.mcmc<-function(x,theta,M)
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


#stepStochasticApproximation
step.SA <-function(x,Z,s.old,gamma)
{
  S<-compute.stat(x,Z)
  s11<-s.old$s1+gamma*(S$s1-s.old$s1)
  s12<-s.old$s2+gamma*(S$s2-s.old$s2)
  s13<-s.old$s3+gamma*(S$s3-s.old$s3)
  s.new<-list(s1=s11,s2=s12,s3=s13)
  return(s.new)
}

step.SA_nesterov <-function(x,Z,s.old,soldd,gamma,sigma)
{ 
  S<-compute.stat(x,Z)
  s11<-s.old$s1+gamma*(S$s1-s.old$s1)+sigma*(s.old$s1-soldd$s1)
  s12<-s.old$s2+gamma*(S$s2-s.old$s2)+sigma*(s.old$s2-soldd$s2)
  s13<-s.old$s3+gamma*(S$s3-s.old$s3)+sigma*(s.old$s3-soldd$s3)
  s.new<-list(s1=s11,s2=s12,s3=s13)
  return(s.new)
}
