library(lme4)
library(rstan)
library(shinystan)#for great model viz
library(ggplot2)#for great viz in general
data(Penicillin)


library(RColorBrewer)
#where the STAN model is saved
#simulate some data
set.seed(20161110)
N<-100 #sample size
J<-10 #number of plant species
id<-rep(1:J,each=10) #index of plant species
K<-3 #number of regression coefficients
#population-level regression coefficient
gamma<-c(2,-1,3)
#standard deviation of the group-level coefficient
tau<-c(0.3,2,1)
#standard deviation of individual observations
sigma<-1
#group-level regression coefficients
beta<-mapply(function(g,t) rnorm(J,g,t),g=gamma,t=tau) 
#the model matrix
X<-model.matrix(~x+y,data=data.frame(x=runif(N,-2,2),y=runif(N,-2,2)))
y<-vector(length = N)
for(n in 1:N){
  #simulate response data
  y[n]<-rnorm(1,X[n,]%*%beta[id[n],],sigma)
}
#run the model
m_hier<-stan(file="pendulum.stan",data=list(N=N,J=J,K=K,id=id,X=X,y=y))

print(m_hier,pars=c("gamma","tau","sigma"))




height = rnorm(100)
earn = height + rnorm(100)

dat = data.frame(height,earn)

qplot(height,earn)+geom_smooth(method=lm,se=F)


# Important annoying fact #1: STAN needs data as a list not a dataframe!!
# specify data as well as meta data (i.e. the number of groups)


earn_dat <- list(N = 100 , #specify number of observations as a scalar
                    earn = earn, # data vector
                    height = height # data vector (predictor) 
                    )

earn_code = 'data {
  // First we declare all of our variables in the data block
  int<lower=0> N;// Number of observations
  vector[N] earn; //Identify our predictor as a vector
  vector[N] height;  //Identify our outcome variable as a vector
}
parameters {
  vector[2] beta; //Our betas are a vector of length 2 (intercept and slope)
  real<lower=0> sigma; //error parameter
}
model {
  //Priors
  beta[1] ~ normal( 5 , .001); //intercept
  beta[2] ~ normal( 0 , 100 ); //slope
  sigma ~ uniform( 0 , 100 ); //error
  earn ~ normal(beta[1] + beta[2] * height, sigma);
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


