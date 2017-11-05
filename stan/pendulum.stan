/*A simple example of an hierarchical model*/
data {
  int<lower=1> N; //the number of observations
  int<lower=1> J; //the number of groups
  int<lower=1> K; //number of columns in the model matrix
  int<lower=1,upper=J> id[N]; //vector of group indeces
  matrix[N,K] X; //the model matrix
  vector[N] y; //the response variable
}
parameters {
  vector[K] gamma; //population-level regression coefficients
  vector<lower=0>[K] tau; //the standard deviation of the regression coefficients

  vector[K] beta[J]; //matrix of group-level regression coefficients
  real<lower=0> sigma; //standard deviation of the individual observations
}
model {
  vector[N] mu; //linear predictor
  //priors
  gamma ~ normal(0,5); //weakly informative priors on the regression coefficients
  tau ~ cauchy(0,2.5); //weakly informative priors, see section 6.9 in STAN user guide
  sigma ~ gamma(2,0.1); //weakly informative priors, see section 6.9 in STAN user guide
  
  for(j in 1:J){
   beta[j] ~ normal(gamma,tau); //fill the matrix of group-level regression coefficients 
  }
  
  for(n in 1:N){
    mu[n] = X[n] * beta[id[n]]; //compute the linear predictor using relevant group-level regression coefficients 
  }

  //likelihood
  y ~ normal(mu,sigma);
}