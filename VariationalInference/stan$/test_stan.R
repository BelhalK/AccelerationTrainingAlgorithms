library(rstan)
library(mvtnorm)

qc1 = 0.0001
deltaT = 0.01
nSamples = 100
m0 = c(1.6, 0)
g = 9.81
t0 = 0.0
ts = seq(deltaT,nSamples * deltaT,deltaT)

bigQ = matrix(c(qc1 * deltaT^3 / 3, qc1 * deltaT^2 / 2,
                qc1 * deltaT^2 / 2,       qc1 * deltaT
                ),
              nrow = 2,
              ncol = 2,
              byrow = TRUE
              )

samples <- stan(file = 'pendulum.stan',
                data = list (
                    T  = nSamples,
                    y0 = m0,
                    t0 = t0,
                    ts = ts,
                    theta = array(g, dim = 1),
                    sigma = c(bigQ[1,1], bigQ[2,2]),
                    refresh = -1
                ),
                algorithm="Fixed_param",
                seed = 42,
                chains = 1,
                iter =1
                )


s <- extract(samples,permuted=FALSE)
plot(s[1,1,1:100])


zStan <- sin(s[1,1,1:nSamples])

estimates <- stan(file = 'penduluminfer.stan',
                  data = list (
                      T  = nSamples,
                      y0 = m0,
                      z  = zStan,
                      t0 = t0,
                      ts = ts
                  ),
                  seed = 42,
                  chains = 1,
                  iter = 1000,
                  warmup = 500,
                  refresh = -1
                  )

e <- extract(estimates,pars=c("theta[1]","sigma[1]","lp__"),permuted=TRUE)

plot(e$lp__[1:100])