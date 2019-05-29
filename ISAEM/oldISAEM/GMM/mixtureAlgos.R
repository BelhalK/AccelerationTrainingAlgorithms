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

mixt.simulate_new <-function(n,weight,mu,sigma)
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

mixt.pdf <- function(df,x,K,n)
{
  col.names <- c("iteration", paste0("Observed data log pdf"))
  pdf <- matrix(NA,K+1,2)
  a <- matrix(NA,n,1)
  for (k in 0:K)
  {
    val <- 0
    for (i in 0:n-1)
    {
      a[i+1] <- -log(sqrt(2*pi)) + log((df[k+1,2]/df[k+1,6])*exp(-(x[i+1]-df[k+1,4])^2/(2*df[k+1,6]^2))+(df[k+1,3]/df[k+1,7])*exp(-(x[i+1]-df[k+1,5])^2/(2*df[k+1,7]^2)))
      val <- sum(val,a[i+1])
    }
    pdf[k+1,] <- c(k, val)
  }

  pdf <- as.data.frame(pdf)
  names(pdf) <- col.names
  return(pdf)
}

mixt.pdftest <- function(df,x,K,n)
{
  col.names <- c("iteration", paste0("Observed data log pdf"))
  pdf <- matrix(NA,K+1,2)
  a <- matrix(NA,n,1)
  for (k in 0:K)
  {
    val <- 0
    for (i in 0:n-1)
    {
      a[i+1] <- log(df$sigma1[k+1]*dnorm(x[i+1],df$mu1[k+1],df$sigma1[k+1])+df$sigma2[k+1]*dnorm(x[i+1],df$mu2[k+1],df$sigma2[k+1]))
      val <- sum(val,a[i+1])
    }
    pdf[k+1,] <- c(k, val)
  }

  pdf <- as.data.frame(pdf)
  names(pdf) <- col.names
  return(pdf)
}

mixt.pdff <- function(df,x,K,n)
{
  col.names <- c("iteration", paste0("Observed data log pdf"))
  pdf <- matrix(NA,K+1,2)
  a <- matrix(NA,n,1)
  for (k in 0:K)
  {
    val <- 0
    for (i in 0:n-1)
    {
      a[i+1] <- log(df$sigma1[k+1]*dnorm(x[i+1],df$mu1[k+1],df$sigma1[k+1])+df$sigma2[k+1]*dnorm(x[i+1],df$mu2[k+1],df$sigma2[k+1])+df$sigma3[k+1]*dnorm(x[i+1],df$mu3[k+1],df$sigma3[k+1]))
      # a[i+1] <- -log(sqrt(2*pi)) + log((df[k+1,2]/df[k+1,6])*exp(-(x[i+1]-df[k+1,4])^2/(2*df[k+1,6]^2))+(df[k+1,3]/df[k+1,7])*exp(-(x[i+1]-df[k+1,5])^2/(2*df[k+1,7]^2))+(df[k+1,4]/df[k+1,10])*exp(-(x[i+1]-df[k+1,7])^2/(2*df[k+1,10]^2)))
      val <- sum(val,a[i+1])
    }
    pdf[k+1,] <- c(k, val)
  }

  pdf <- as.data.frame(pdf)
  names(pdf) <- col.names
  return(pdf)
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

#-------------------------------------
mixt.em_replace <- function(x, theta0, K, nb_r)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  tau<-compute.tau(x,theta)
  for (k in 1:K)
  {
    i<-sample(1:n,nb_r)
    tau[i,] <- compute.tau_replace(x[i],theta,nb_r)
    s <- compute.stat(x,tau)
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

#-------------------------------------
mixt.mcem1 <- function(x, theta0, K, K1=1, M=c(1,10))
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  K2 <- K - K1 
  MK <-c(rep(M[1],K1),rep(M[2],K2))

  # print(c(MK[K],sum(MK)))
  
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  for (k in 1:K)
  {
    Z<-step.S(x,theta,MK[k])
    s<-step.SA(x,Z,s,1)
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}



mixt.mcem2 <- function(x, theta0, K, K1=NULL, M=1, alpha=1)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  if (is.null(K1))  K1 <- 1
  K2 <- K - K1 

  if (length(M)==1)
  MK <-c(rep(M,K1),M+round((1:K2)^alpha))
  else
  MK <-c(rep(M[1],K1),rep(M[2],K2))
  # print(c(MK[K],sum(MK)))
  
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  for (k in 1:K)
  {
    Z<-step.S(x,theta,MK[k])
    s<-step.SA(x,Z,s,1)
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

mixt.mcem3 <- function(x, theta0, K, K1=1, M=c(10,1), alpha=1)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  K2 <- K - K1 
  M <-c(rep(M[1],K1),rep(M[2],K2))
  # print(c(MK[K],sum(MK)))
  
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  gamma<-c(rep(1,K1),1/(1:K2)^alpha)
  
  for (k in 1:K)
  {
    Z<-step.S(x,theta,M[k])
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

mixt.mcem4 <- function(x, theta0, K, K1=1, M=c(10,1), alpha=0.5)
{
  G<-length(mu)
  col.names <- c(paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  K2 <- K - K1
  gamma<-c(rep(1,K1),1/(1:K2)^alpha)
  delta<-c(rep(1,K1),1/(1:K2))
  M <-c(rep(M[1],K1),rep(M[2],K2))
  
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  M0 <- matrix(rep(0,G),nrow=1)
  # s.mean <- list(s1=M0, s2=M0, s3=M0)
  t.mean <- list(p=M0, mu=M0, sigma=M0)
  for (k in 1:K)
  {
    #    M <- ceiling(k/5)
    Z<-step.S(x,theta,M[k])
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    km <- delta[k]
    for (j in (1:3))
    {
      # s.mean[[j]] <- rbind(s.mean[[j]],s[[j]]*delta[k] + s.mean[[j]][k,]*(1-delta[k]) )
      t.mean[[j]] <- rbind(t.mean[[j]],theta[[j]]*delta[k] + t.mean[[j]][k,]*(1-delta[k]) )
    }
  }
  # theta<-step.M(s.mean,n)
  theta<-t.mean
  df <- as.data.frame(theta)
  df[1,] <- c(theta0$p, theta0$mu, theta0$sigma)
  names(df) <- col.names
  df <- cbind(data.frame(iteration = 0:K),df)
  return(df)
}


mixt.saem1 <- function(x, theta0, K, K1=NULL, M=1, alpha=1)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  if (is.null(K1))  K1 <- 1
  K2 <- K - K1 #second phase iterations
  if (length(alpha)==1)
  gamma<-c(rep(1,K1),1/(1:K2)^alpha)
  else{
    L <- 10
    KL <- round(K2/L)
    alpha <- seq(alpha[1], alpha[2], length.out = L)
    gamma <- rep(1,K1)
    dl <- 0
    for (l in (1:L))
    {
      if (l==L)  KL <- K2 - (L-1)*KL
      gamma <- c(gamma,1/(dl + (1:KL))^alpha[l])
      dl <- (dl + KL)^(alpha[l]/alpha[l+1])
    }
  }
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  
  
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  for (k in 1:K)
  {
    # if (k %% 1000==0)
    # {
    #   print(k)
    # }
    Z<-step.S(x,theta,M)
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

mixt.saem1_nesterov <- function(x, theta0, K, K1=NULL, M=1, alpha=1, sigma=0.6)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  if (is.null(K1))  K1 <- 1
  K2 <- K - K1 #second phase iterations
  if (length(alpha)==1)
  gamma<-c(rep(1,K1),1/(1:K2)^alpha)
  else{
    L <- 10
    KL <- round(K2/L)
    alpha <- seq(alpha[1], alpha[2], length.out = L)
    gamma <- rep(1,K1)
    dl <- 0
    for (l in (1:L))
    {
      if (l==L)  KL <- K2 - (L-1)*KL
      gamma <- c(gamma,1/(dl + (1:KL))^alpha[l])
      dl <- (dl + KL)^(alpha[l]/alpha[l+1])
    }
  }
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  
  
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  sold<-list(s1=0,s2=0,s3=0)
  for (k in 2:K)
  {
    # if (k %% 1000==0)
    # {
    #   print(k)
    # }
    Z<-step.S(x,theta,M)
    sold<-step.SA(x,Z,s,gamma[k-1])
    s<-step.SA_nesterov(x,Z,s,sold,gamma[k],sigma)
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

mixt.saem1_replace1 <- function(x, theta0, K, K1=NULL, M=1, alpha=1, nb_r=1)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  if (is.null(K1))  K1 <- 1
  K2 <- K - K1 #second phase iterations
  if (length(alpha)==1)
  gamma<-c(rep(1,K1),1/(1:K2)^alpha)
  else{
    L <- 10
    KL <- round(K2/L)
    alpha <- seq(alpha[1], alpha[2], length.out = L)
    gamma <- rep(1,K1)
    dl <- 0
    for (l in (1:L))
    {
      if (l==L)  KL <- K2 - (L-1)*KL
      gamma <- c(gamma,1/(dl + (1:KL))^alpha[l])
      dl <- (dl + KL)^(alpha[l]/alpha[l+1])
    }
  }

  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  Z<-step.S_replace(x,theta,M)
  n<-length(x)
  for (k in 1:K)
  {
    # if (k %% 1000==0)
    # {
    #   print(k)
    # }
    i<-sample(1:n,nb_r)
    Z[i,,]<-step.S_replace(x[i],theta,M)
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.saem1_replace_K2 <- function(x, theta0, K, K1=NULL, M=1, alpha=1, nb_r=1)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  if (is.null(K1))  K1 <- 1
  K2 <- K - K1 #second phase iterations
  if (length(alpha)==1)
  gamma<-c(rep(1,K1),1/(1:K2)^alpha)
  else{
    L <- 10
    KL <- round(K2/L)
    alpha <- seq(alpha[1], alpha[2], length.out = L)
    gamma <- rep(1,K1)
    dl <- 0
    for (l in (1:L))
    {
      if (l==L)  KL <- K2 - (L-1)*KL
      gamma <- c(gamma,1/(dl + (1:KL))^alpha[l])
      dl <- (dl + KL)^(alpha[l]/alpha[l+1])
    }
  }
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  Z<-step.S_replace(x,theta,M)
  n<-length(x)
  for (k in 1:K)
  { 
    if(k<K1){
    Z<-step.S(x,theta,M)
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
    }
    else{
    i<-sample(1:n,nb_r)
    Z[i,,]<-step.S_replace(x[i],theta,M)
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
    }
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.saem2 <- function(x, theta0, K, K1=NULL, M=1, alpha=0.5)
{
  G<-length(mu)
  col.names <- c(paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  if (is.null(K1))  K1 <- 1
  K2 <- K - K1
  gamma<-c(rep(1,K1),1/(1:K2)^alpha)
  delta<-c(rep(1,K1),1/(1:K2))
  
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  M0 <- matrix(rep(0,G),nrow=1)
  s.mean <- list(s1=M0, s2=M0, s3=M0)
  for (k in 1:K)
  {
    #    M <- ceiling(k/5)
    Z<-step.S(x,theta,M)
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    km <- delta[k]
    for (j in (1:3))
    {
      s.mean[[j]] <- rbind(s.mean[[j]],s[[j]]*delta[k] + s.mean[[j]][k,]*(1-delta[k]) )
    }
  }
  theta<-step.M(s.mean,n)
  df <- as.data.frame(theta)
  df[1,] <- c(theta0$p, theta0$mu, theta0$sigma)
  names(df) <- col.names
  df <- cbind(data.frame(iteration = 0:K),df)
  return(df)
}

mixt.saem3 <- function(x, theta0, K, K1=NULL, M=c(10,1), alpha=0.5)
{
  G<-length(mu)
  col.names <- c(paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  if (is.null(K1))  K1 <- 1
  K2 <- K - K1
  gamma<-c(rep(1,K1),1/(1:K2)^alpha)
  delta<-c(rep(1,K1),1/(1:K2))
  M <-c(rep(M[1],K1),rep(M[2],K2))
  
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  M0 <- matrix(rep(0,G),nrow=1)
  t.mean <- list(p=M0, mu=M0, sigma=M0)
  for (k in 1:K)
  {
    #    M <- ceiling(k/5)
    Z<-step.S(x,theta,M[k])
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    km <- delta[k]
    for (j in (1:3))
    {
      t.mean[[j]] <- rbind(t.mean[[j]],theta[[j]]*delta[k] + t.mean[[j]][k,]*(1-delta[k]) )
    }
  }
  theta<-t.mean
  df <- as.data.frame(theta)
  df[1,] <- c(theta0$p, theta0$mu, theta0$sigma)
  names(df) <- col.names
  df <- cbind(data.frame(iteration = 0:K),df)
  return(df)
}


mixt.saem1_mcmc <- function(x, theta0, K, K1=NULL, M=1, alpha=1)
{
  G<-length(mu)
  col.names <- c("iteration", paste0("p",1:G), paste0("mu",1:G), paste0("sigma",1:G))
  
  if (is.null(K1))  K1 <- 1
  K2 <- K - K1 #second phase iterations
  if (length(alpha)==1)
  gamma<-c(rep(1,K1),1/(1:K2)^alpha)
  else{
    L <- 10
    KL <- round(K2/L)
    alpha <- seq(alpha[1], alpha[2], length.out = L)
    gamma <- rep(1,K1)
    dl <- 0
    for (l in (1:L))
    {
      if (l==L)  KL <- K2 - (L-1)*KL
      gamma <- c(gamma,1/(dl + (1:KL))^alpha[l])
      dl <- (dl + KL)^(alpha[l]/alpha[l+1])
    }
  }
  theta.est <- matrix(NA,K+1,3*G+1)
  theta.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  
  
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  for (k in 1:K)
  {
    # if (k %% 1000==0)
    # {
    #   print(k)
    # }
    Z<-step.mcmc(x,theta,M)
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}