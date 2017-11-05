library(lme4)
library(rstan)
library(shinystan)#for great model viz
library(ggplot2)#for great viz in general
data(Penicillin)

library(RColorBrewer)

age = rnorm(100)

base_0 <- 1
slope_0 <- 1
height = base_0 +slope_0*age + rnorm(100)


dat = data.frame(age,height)

qplot(age,height)+geom_smooth(method=lm,se=F)


# Important annoying fact #1: STAN needs data as a list not a dataframe!!
# specify data as well as meta data (i.e. the number of groups)


earn_dat <- list(N = 100 , #specify number of observations as a scalar
                    height = height, # data vector
                    age = age # data vector (predictor) 
                    )

earn_code = 'data {
  // First we declare all of our variables in the data block
  int<lower=0> N;// Number of observations
  vector[N] age; //Identify our predictor as a vector
  vector[N] height;  //Identify our outcome variable as a vector
}
parameters {
  vector[2] beta; //Our betas are a vector of length 2 (intercept and slope)
  real<lower=0> sigma; //error parameter
}
model {
  //Priors
  beta[1] ~ normal( 1 , .001); //intercept
  beta[2] ~ normal( 0 , 0.1 ); //slope
  sigma ~ uniform( 0 , 1 ); //error
  height ~ normal(beta[1] + beta[2] * age, sigma);
}'


fit1 <- stan(model_code = earn_code, data = earn_dat,
             warmup = 100,
             iter = 1000, 
             chains = 4)

print(fit1)


# extract posterior samples for each parameter
fit1_samples = extract(fit1)
str(fit1_samples)
# subset just the betas
betas = fit1_samples[[1]]

qplot(betas[,1]) # intercept posterior samples
qplot(betas[,2]) # slope posterior samples


