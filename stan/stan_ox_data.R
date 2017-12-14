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


model <- 'data {
          int<lower=0> N;// Number of observations
          real beta1_pop;
          real beta2_pop;
          real<lower=0>  pres;
          real<lower=0> omega_beta1;
          real<lower=0> omega_beta2;
          vector[N] age; //predictor
          vector[N] height;  //response
        }
        parameters {
          vector[2] beta;
        }
        model {
          //Priors
          beta[1] ~ normal( beta1_pop , omega_beta1);
          beta[2] ~ normal( beta2_pop , omega_beta2);
          for (n in 1:N)
            height[n] ~ normal(beta[1] + beta[2] * age[n], pres);
        }'



modelstan <- stan_model(model_name = "oxboys",model_code = model)

# fit <- stan(model_code = earn_code, data = earn_dat,
#              warmup = 1,
#              iter = 3, 
#              chains = 1)

stan_data <- list(N = length(Dargs$yobs[Dargs$IdM==i]),height = Dargs$yobs[Dargs$IdM==i]
                ,age = Dargs$XM[Dargs$IdM==i,],
                beta1_pop=mean.phiM[i,1],beta2_pop=mean.phiM[i,2],
                omega_beta1=omega.eta[1,1],omega_beta2=omega.eta[2,2],
                pres=1)

fit <- sampling(stan.model, data = stan_data,algorithm = "NUTS", chains = 1,iter = 100, warmup = 1)
print(fit)

# extract posterior samples for each parameter
fit_samples = extract(fit)
str(fit_samples)
# subset just the betas
betas = fit_samples[[1]]

qplot(betas[,1]) # intercept posterior samples
qplot(betas[,2]) # slope posterior samples


