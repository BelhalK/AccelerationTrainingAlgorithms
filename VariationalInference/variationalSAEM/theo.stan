data {
int<lower=0> N;
vector[N] y;
vector[N] time;
}
parameters {
real alpha;
real alpha_pop;
real beta;
real beta_pop;
real omega;
real<lower=0> sigma;
}
model {
alpha ~ normal(alpha_pop,omega);
beta ~ normal(beta_pop,omega);
y ~ normal(alpha + beta * time, sigma);
}