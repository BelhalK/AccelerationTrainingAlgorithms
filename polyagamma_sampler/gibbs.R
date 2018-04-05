# Here will be some functions for Polya-Gamma Gibbs sampling

# Naive sampler for PG(1, z)
# (based on the finite approximation of infinite sum)
rpolyagamma_naive = function(z, n_terms = 100){
  g = rexp(n_terms, 1)
  out = 1 / (2*pi**2) * sum(g / ((1:n_terms - 1/2)**2 + z**2 / (4*pi**2)))
  return(out)
}#end function


#FUNCTION: cdf(x) of inverse Normal distribution with parameters mu and lambda
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
pinversen = function(x,mu,lambda){
  #check if x, mu and lambda are positive real numbers
  if ((x<=0)|(mu<=0)|(lambda<=0)){
    stop("Parameters in pinversen() are not of the correct type");
  }#end if

  #assign variables
  sqrt_lambda_over_x = sqrt(lambda/x);
  x_over_mu = x/mu;

  #work out the cdf and return it
  cdf = pnorm(sqrt_lambda_over_x*(x_over_mu-1));
  if (cdf!=1){
    cdf = cdf + exp(2*lambda/mu)*pnorm(-sqrt_lambda_over_x*(x_over_mu+1));
  }#end if
  return (cdf);
}#end pinversen


#FUNCTION: Piecewise coefficients
#PARAMETERS:
#Parameter x, nth term, truncate at t
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
a_n = function(x,n,t){
  #check if n is positive integer, x is positive real, t is positive real
  if (any(x<=0)|(t<=0)|(round(n)!=n)|(n<0)){
    stop("Parameters in a_n are not of the correct type");
  }#end if

  #for a <= t
  a = (x<=t)*(pi*(n+0.5)*(sqrt(2/(pi*x)))^3*exp(-2*(n+0.5)^2/x));
  #for a > t
  a = a + (x>t)*(pi*(n+0.5)*exp(-0.5*(n+0.5)^2*pi^2*x));

  #return a
  return(a);
}#end a_n

#FUNCTION: Generate sample from inverse Normal with parameters (mu, 1) truncated with a max of t
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
rinversen = function(mu,t){
  #check if mu and t are positive real
  if ((mu<=0)|(t<=0)){
    stop("Parameters in rinversen are not of the correct type");
  }#end if

  #x is the random sample
  x = t;

  #while x is equal or more than t, sample
  while(x>=t){

    #for large mu, use chi-squared approximation
    if (mu>t){
      repeat{
        #sample exponential (full details in paper)
        repeat{
          E = rexp(2);
          if (E[1]^2<=2*E[2]/t){
            E = E[1];
            break;
          }#end if
        }#end repeat

        #assign x, alpha
        x = t/(1+t*E)^2;
        alpha = exp(-0.5*x/mu^2);

        #accept if U(0,1) <= alpha
        if (runif(1)<=alpha){
          break;
        }#end if

      }#end repeat
    }#end if

    #else for small mu...
    else{
      y = (rnorm(1))^2;
      x = mu+0.5*mu*mu*y-0.5*mu*sqrt(4*mu*y+(mu*y)^2);
      if (runif(1)>mu/(mu+x)){
        x = mu*mu/x;
      }#end if
    }#else

  }#end while

  #return the sample
  return(x);
}#end rinverse

#DEBUG FUNCTION: pdf(x) of inverse normal with parameter mu, lambda
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
dinversen = function(x,mu,lambda){
  return(sqrt(lambda/(2*pi*x^3))*exp(-lambda*(x-mu)^2/(2*mu^2*x)));
}#end dinversen

#DEBUG FUNCTION: pdf(x|x<t) where x~IG(mu,lambda)
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
dtruncatedinversen = function(x,mu,lambda,t){
  return(dinversen(x,mu,lambda)/pinversen(t,mu,lambda));
}#end dtruncatedinversen

#DEBUG FUNCTION: plot histogram and pdf of truncated IG(mu,1) at t
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
testinversen = function(mu,t){
  x = replicate(10000,rinversen(mu,t));
  #Observeable: the pdf can't be evaluated for mu<0.0029
  x_plot = seq(from=min(x),to=max(x),length=100000);
  f_plot = sapply(x_plot,dtruncatedinversen,mu=mu,lambda=1,t=t);
  hist(x,freq=FALSE);
  lines(x_plot,f_plot);
}#end testinversen


#FUNCTION: Generate sample from PolyaGamma(a,z) with truncation t
#PARAMETERS: a is an integer, z and t are positive real numbers
#AUTHOR: SHERMAN IP
#DATE: 16/10/15 - 19/10/15
#NOTES: Fails for high z
rpolyagamma = function(a,z,t){
  #check if a is integer, t are positive real numbers
  if ((a!=round(a))|(a<=0)|(t<=0)){
    stop("Parameters in rpolyagamma are not of the correct type");
  }#end if

  #if a is not 1, ie >1, then add together rpolyagamma sample a times
  if (a!=1){
    x= 0; #assign x
    #a times
    for (i in seq(from=1,to=a)){
      #add sample from polyagamma(1,z)
      x = x + rpolyagamma(1,z,t);
    }#end for
    #return x
    return (x);
  }#end if

  #else a is 1
  else{
    z = abs(z)/2;
    mu = 1/z;
    K = pi*pi/8+z*z/2;

    #calculate mixture coefficient
    p = pi/(2*K)*exp(-K*t);
    q = 2*exp(-z)*pinversen(t,mu,lambda=1);
    r = p/(p+q);
    #if the mixture coefficent is not finite, stop the program
    if (!is.finite(r)){
      stop("Mixture coefficients are not finite in rpolyagamma");
    }#end if

    #put in global REJECTION count
    REJECTION <<- 0;

    #accept-reject sample, repeat until accept
    repeat{
      #sample x from mixture model

      #probability p/(p+q), sample truncated exp
      if(runif(1)<r){
        x = t+rexp(1)/K;
      }#end if
      #else sample from inverse n
      else{
        x = rinversen(mu,t);
      }#end else

      #get 0th coefficient
      S = a_n(x,0,t);
      y = runif(1)*S;
      n = 0;
      repeat{
        n = n+1;
        #for odd n
        if ((n%%2)==1){
          S = S-a_n(x,n,t);
          #if y is smaller than S, accept and return it
          if (y<S){
            return(x/4);
          }#end if
        }#end if
        #for even n
        else{
          S = S+a_n(x,n,t);
          #if y is bigger than S, reject it and increase the global counter
          if (y>S){
            REJECTION <<- REJECTION + 1;
            break;
          }#end if
        }#end else
      }#end repeat
    }#end repeat
  }#end else
}#end rpolyagamma

#PLOT FUNCTION
#plot the histogram and curve of samples and pdf of polya gamma used for the report
plotHistogramPolyaGamma = function(){
  par(mfrow=c(2,3));
  testpolyagamma(0.1,0.64);
  testpolyagamma(1,0.64);
  testpolyagamma(10,0.64);
  testpolyagammanaive(10,2);
  testpolyagammanaive(10,10);
  testpolyagammanaive(10,100);
}#end plotHistogramPolyaGamma

#DEBUG FUNCTION: plot histogram and pdf of PG(1,z) with truncation t
#AUTHOR: SHERMAN IP
#DATE: 16/10/15
testpolyagamma = function(z,t){
  x = replicate(10000,rpolyagamma(1,z,t));

  #Observeable: the pdf can't be evaluated for mu<0.0029
  x_plot = seq(from=min(x),to=max(x),length=1000);
  f_plot = sapply(x_plot,dpolyagamma,z=z,t=t);
  hist(x,freq=FALSE,main=paste("Proposed PG(1,",toString(z),") sampler",sep=""),xlim=c(0,max(x)), ylim=c(0, max(f_plot)*1.2),breaks=20);
  lines(x_plot,f_plot);
}#end testinversen


#DEBUG FUNCTION: plot histogram and pdf of PG(1,z) using the naive sampler
#AUTHOR: SHERMAN IP
#DATE 20/10/15
testpolyagammanaive = function(z,n){
  x = replicate(10000,rpolyagamma_naive(z,n_terms = n));
  #Observeable: the pdf can't be evaluated for mu<0.0029
  x_plot = seq(from=min(x),to=max(x),length=1000);
  f_plot = sapply(x_plot,dpolyagamma,z=z,t=t);
  hist(x,freq=FALSE,main=paste("Naive PG(1,",toString(z),") sampler with ",toString(n)," term(s)",sep=""),xlim=c(0,max(x)),ylim=c(0, max(f_plot)*1.2),breaks=20);
  lines(x_plot,f_plot);
}#end testinversen


#DEBUG FUNCTION: pdf of PG(1,z) at x and truncation point t
#AUTHOR: SHERMAN IP
#DATE: 16/10/15
dpolyagamma = function(x,z,t){
  pdf = 0;
  #use the first 100 terms in infinite sum
  for (i in 0:100){
    pdf = pdf + (-1)^i*a_n(4*x,i,t);
  }#end for
  return(4*cosh(z/2)*exp(-z*z*4*x/8)*pdf);
}#end dpolyagamma

#EXPERIMENT FUNCTION: investegate the rejection rate of the polyagamma sampler
#PARAMETERS: sample PG(1,10^b_exp_array) n times
rejectionPolyaGamma = function(b_exp_array,n){
  #check if b is a vector and n is a positive integer
  if ( (!is.vector(b_exp_array)) | (n!=round(n)) | (n<=0) ){
    stop("Paramters in rejectionPolyaGamma are not of the correct type");
  }
  rejection_matrix = matrix(0,nrow=n,ncol=length(b_exp_array));
  for (j in seq_len(length(b_exp_array))){
    b = 10^(b_exp_array[j]);
    for (i in 1:n){
      rpolyagamma(1,b,0.64);
      rejection_matrix[i,j] = REJECTION;
    }#end for
  }#end for

  #plot the training and testing error
  boxplot(rejection_matrix,names=paste("10E",sapply(b_exp_array,toString),sep=""),xlab="b",ylab="Number of rejection");
  mean = colMeans(rejection_matrix);
  errbar(b_exp_array,mean,apply(rejection_matrix,2,min),apply(rejection_matrix,2,max));

}#end rejectionPolyaGamma

# Generate parameter vector beta from a multivariate normal
# beta ~ N(m, V), see details in the paper
generate_mv_normal = function(w, y, X, b, Binv){
  temp = t(X) %*% (w * X) + Binv
  kappa = y - 0.5
  V = chol2inv(chol(temp))
  m = V %*% as.vector(t(X) %*% kappa + Binv %*% b)
  # now generate beta ~ mvNorm(m, V)
  beta = m + t(chol(V)) %*% rnorm(length(b))
  return(beta)
}

# Check whether the user has provided arguments
# with matching dimensionalities
check_dimensions = function(y, X, b, B){
  if(!is.vector(y)) stop("y must be a vector")
  if(!is.vector(b)) stop("b must be a vector")
  if(!is.matrix(X)) stop("X must be a matrix")
  if(!is.matrix(B)) stop("B must be a matrix")
  if(length(y) != nrow(X)) stop("nrow(X) must equal length(y)")
  if(length(b) != ncol(X)) stop("ncol(X) must equal length(b)")
}

#' Gibbs two-step sampling procedure
#' for parameter vector beta and the latent Polya-Gamma variables
#'
#' @param y binary vector of observations
#' @param X design matrix of covariates
#' @param lambda The diagonal element of the precision matrix of the prior distribution (see also parameter B)
#' @param b The prior mean of the prior distribution. Defaults to a vector of zeros.
#' @param B The prior precision of the prior distribution. Defaults to lambda * identity.
#' @param n_iter The total number of iterations in the MCMC chain
#' @param naive Should the naive approximation be used to generate the Polya-Gamma distribution
#' @param naive_n_terms If the naive approximation is used, then this specifies number of terms in the finite sum.
#' @param t The parameter in the accept-reject algorithm for sampling from Polya-Gamma distribution (see paper for details).
#'
#' @return list containing the MCMC samples from the posterior distribution of beta
#' @examples
#' data = generate_from_simple_logistic_model(n=100)
#' obj = gibbs_sampler(data$y, data$X, lambda=0.001, n_iter_total=100, burn_in=50)
#' plot(obj)
#' @export gibbs_sampler
#'
gibbs_sampler = function(y, X, lambda = 0.0001, b=rep(0, ncol(X)), B=lambda*diag(ncol(X)), n_iter_total = 200, burn_in = 100, naive = FALSE, naive_n_terms = 100, t = 0.64){
  # Check if everything is OK with dimensions
  check_dimensions(y, X, b, B)
  # number of parameters
  m = ncol(X)
  # number of data points
  n = nrow(X)
  # Starting values for beta; initialise w
  beta = b
  w = rep(NA, n)

  # Store the values of all betas and all w
  beta_all = matrix(0, n_iter_total, m)
  w_all = matrix(0, n_iter_total, n)

  # initialise progressbar
  pb = txtProgressBar(min = 0, max = n_iter_total, initial = 0)

  for(k in 1:n_iter_total){
    # draw elements of w from PG
    for(i in 1:n){
      psi = as.numeric(X[i, ] %*% beta)
      if(naive) w[i] = rpolyagamma_naive(psi, naive_n_terms) else w[i] = rpolyagamma(1, psi, t)
    }
  
    # draw beta from a multivariate normal
    beta = generate_mv_normal(w, y, X, b, B)
    beta_all[k, ] = beta
    w_all[k, ] = w
    if(k%%50==0) setTxtProgressBar(pb, k)
  }
  close(pb)
  selected_ind = (burn_in+1):n_iter_total
  out = list("beta" = beta_all[selected_ind, ], "w" = w_all[selected_ind, ], "burn_in" = burn_in)
  class(out) = "PG"
  return(out)
}

#' Print a summary of the MCMC chains
#'
#' @export
print.PG = function(obj){
  beta = obj$beta
  posterior_mean = round(colMeans(beta), 3)
  posterior_sd = round(apply(beta, 2, sd), 3)
  s = "
  MCMC sample from the posterior distribution of beta.
  Chain length: %d (with %d burn-in removed).
  Number of parameters: %d.
  Posterior means: %s.
  Posterior standard deviations: %s.
  "
  cat(sprintf(s,
              nrow(beta),
              obj$burn_in,
              ncol(beta),
              paste(posterior_mean, collapse=", "),
              paste(posterior_sd, collapse=", ")))
}


#' Traceplot and autocorrelation plot for the beta chain
#' @param which_parameters Vector of indexes which components of beta vector should be plotted.
#'
#' @export
plot.PG = function(obj, which_parameters = 1:ncol(obj$beta)){
  X = obj$beta[, which_parameters, drop=FALSE]

  layout_mat = create_layout_matrix(ncol(X))
  layout(layout_mat)
  for(j in 1:ncol(X)){
    x = X[, j]
    # traceplot
    plot(x, type="l")
    # autocorrelation plot
    acf(x, main="")
  }
  # restore layout to default
  layout(1)
}

create_layout_matrix = function(i){
  col1 = seq(1, 2*i, 2)
  col2 = col1
  col3 = col1 + 1
  return(cbind(col1, col2, col3))
}