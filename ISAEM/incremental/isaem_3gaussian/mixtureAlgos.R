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
    tau<-step.E(x,theta)
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  
  df <- as.data.frame(theta.est)
  names(df) <- col.names
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
    # if (k%%(n/nb_r) == 1)
    # {
    #   i<-1:nb_r
    # }
    i<-sample(1:n,nb_r)
    Z[i,,]<-step.S_replace(x[i],theta,M)
    s<-step.SA(x,Z,s,gamma[k])
    # i <- i+nb_r
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}




mixt.isaem1 <- function(x, theta0, K, K1=NULL, M=1, alpha=1, nb_r=1,target)
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
  thetat.est <- matrix(NA,K+1,3*G+1)
  thetat.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  st<-list(s1=0,s2=0,s3=0)
  Z<-step.S_replace(x,theta,M)
  Z_old<-Z
  n<-length(x)
  for (k in 1:K)
  {
    diff <- matrix(NA,n,2)
    for (i in 1:n)
    {
      Z_old[i,,] <- Z[i,,]
      Z[i,,]<-step.S_replace(x[i],theta,M)
      st<-step.SA(x,Z,st,gamma[k])
      thetat<-step.M(st,n)
      thetat.est[k+1,] <- c(k, thetat$p, thetat$mu, thetat$sigma)
      dff <- as.data.frame(thetat.est)
      names(dff) <- col.names
      diff[i,2] <- norm(t((dff[(k+1),2:(3*G+1)] - target)))
      diff[i,1] <- i
      Z[i,,] <- Z_old[i,,]
    }
    diff <- diff[order(diff[,2],decreasing=FALSE),]
    i <- diff[1:nb_r,1]
    Z[i,,]<-step.S_replace(x[i],theta,M)
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}

mixt.isaem2 <- function(x, theta0, K, K1=NULL, M=1, alpha=1, nb_r=1,target)
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
  thetat.est <- matrix(NA,K+1,3*G+1)
  thetat.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  st<-list(s1=0,s2=0,s3=0)
  # Z<-step.S_replace(x,theta,M)
  obj <- step.S_isaem2(x,theta,M)
  Z<-obj$var
  tau<-obj$prob
  tau_old <- tau
  Z_old<-Z
  n<-length(x)
  for (k in 1:K)
  {
      # browser()
    diff <- matrix(NA,n,2)
    obj <- step.S_isaem2(x,theta,M)
    tau_test <- obj$prob
    # diff[,2] <- norm(t(tau_test[,] - tau_old[,]))

    if (length(theta0$p) == 2) {
       diff[,2] <- (tau_test[,1] - tau_old[,1])^2 + (tau_test[,2] - tau_old[,2])^2
    } else {
       diff[,2] <- (tau_test[,1] - tau_old[,1])^2 + (tau_test[,2] - tau_old[,2])^2 + (tau_test[,3] - tau_old[,3])^2
    }
    
    diff[,1] <- 1:n
    diff <- diff[order(diff[,2],decreasing=TRUE),]

    i <- diff[1:nb_r,1]
    Z[i,,] <- obj$var[i,,]
    # Z[i,,]<-step.S_isaem2(x[i],theta,M)$var
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
    tau_old <- tau
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}


mixt.isaem3 <- function(x, theta0, K, K1=NULL, M=1, alpha=1, nb_r=1)
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
  thetat.est <- matrix(NA,K+1,3*G+1)
  thetat.est[1,] <- c(0, theta0$p, theta0$mu, theta0$sigma)
  theta<-theta0
  s<-list(s1=0,s2=0,s3=0)
  st<-list(s1=0,s2=0,s3=0)
  liste <- vector("list", length=nsim)
  Z<-step.S_replace(x,theta,M)
  Z_old<-Z
  n<-length(x)
  grade <- matrix(NA,n,2)
  theta1 <- matrix(NA,n,3*G)
  df1 <- as.data.frame(theta1)
  for (k in 1:K)
  {
    diff <- matrix(NA,n,2)
    if (k<2)
    {
      for (l in 1:n)
      {
        Z_old[l,,] <- Z[l,,]
        Z[l,,]<-step.S_replace(x[l],theta,M)
        st<-step.SA(x,Z,st,gamma[k])
        thetat<-step.M(st,n)
        thetat.est[k+1,] <- c(k, thetat$p, thetat$mu, thetat$sigma)
        df <- as.data.frame(thetat.est)
        names(df) <- col.names
        df1[l,] <- df[(k+1),(3*G+1)]
        grade[l,2] <- norm(t((df1[l,] - df[k,(3*G+1)])))
        grade[l,1] <- l
        Z[l,,] <- Z_old[l,,]
      }
      grade <- grade[order(grade[,2],decreasing=TRUE),]
      i <- grade[1:nb_r,1]
      liste[[k]] <- i
      Z[i,,]<-step.S_replace(x[i],theta,M)
      s<-step.SA(x,Z,s,gamma[k])
      theta<-step.M(s,n)
      df <- as.data.frame(theta.est)
      names(df) <- col.names
      theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
    } else
    {
      theta_old <- theta  
      grade <- grade[order(grade[,2],decreasing=TRUE),]
      i <- grade[1:nb_r,1]
      liste[[k]] <- i
      # i <- sample(1:n,nb_r)
      Z[i,,]<-step.S_replace(x[i],theta,M)
      s<-step.SA(x,Z,s,gamma[k])
      theta<-step.M(s,n)
      theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
      df <- as.data.frame(theta.est)
      names(df) <- col.names
      grade[i,2] <- norm(t((df[(k+1),(3*G+1)] - df[k,(3*G+1)])))
      #computing the parameter as if we had chosen the individual liste[k-1]
      Z_old[liste[[k-1]],,]<-step.S_replace(x[liste[[k-1]]],theta_old,M)
      s<-step.SA(x,Z_old,s,gamma[k])
      thetat<-step.M(s,n)
      thetat.est[k+1,] <- c(k, thetat$p, thetat$mu, thetat$sigma)
      dft <- as.data.frame(thetat.est)
      grade[liste[[k-1]],2] <- norm(t((dft[(k+1),(3*G+1)] - df[(k+1),(3*G+1)])))
    }
  }
  # df <- as.data.frame(theta.est)
  # names(df) <- col.names
  return(df)
}





# invcdf <- function(u, m) {
#     return(sqrt(m^2/(1 - (1 - m^2) * u)))
# }



mixt.isaem_bouchard <- function(x, theta0, K, K1=NULL, M=1, alpha=1, nb_r=1)
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
  tau<-step.S_isaem2(x,theta,M)$prob
  n<-length(x)
  m<- 0.2

  for (k in 1:K)
  {
    

    F <- function(x,m) 1-(1/m+x)*exp(-m*x)
    F.inv <- function(y,m){uniroot(function(x){F(x,m)-y},interval=c(0,100))$root}
    F.inv <- Vectorize(F.inv)


    X <- runif(1000,0,1)   # random sample from U[0,1]
    D <- F.inv(X,m)
    # sample1 <- sapply(runif(100), invcdf, m = m)
    a <- hist(D, plot=FALSE, breaks=c(seq(0,40,length=30),Inf))
    l <- rmultinom(n=nb_r,size=99,prob=c(max(a$density),1-max(a$density)))
    # l <- rmultinom(n=nb_r,size=99,prob=c(sample1[1],1-sample1[1]))
    i <- unique(l[1,1:nb_r])
    
    # Z[l[1],,]<-step.S_replace(x[l[1]],theta,M)
    Z[i,,]<-step.S_replace(x[i],theta,M)
  
    m <- m + gamma[k]*tau[41,2]
    # print(m)
    s<-step.SA(x,Z,s,gamma[k])
    theta<-step.M(s,n)
    theta.est[k+1,] <- c(k, theta$p, theta$mu, theta$sigma)
  }
  df <- as.data.frame(theta.est)
  names(df) <- col.names
  return(df)
}