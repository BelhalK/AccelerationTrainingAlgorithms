library(sgd)
library(ggplot2)

generate.data <- function(N, d) {
  l2 <- function(x) sqrt(sum(x**2))
  X <- matrix(rnorm(N*d, mean=0, sd=1/sqrt(N)), nrow=N, ncol=d)
  theta <- runif(d)
  theta <- theta * 6 *sqrt(d) / l2(theta)

  # noise
  ind <- rbinom(N, size=1, prob=.95)
  epsilon <- ind * rnorm(N) + (1-ind) * rep(10 ,N)

  Y <- X %*% theta + epsilon
  return(list(y=Y, X=X, theta=theta))
}

# Dimensions
N <- 1000
d <- 200

# Generate data.
set.seed(42)
data <- generate.data(N, d)
dat <- data.frame(y=data$y, x=data$X)

sgd.theta <- sgd(y ~ .-1, data=dat, model="m", sgd.control=list(method="sgd",
  lr.control=c(15, NA, NA, 1/2), npass=50, pass=T))

plot(sgd.theta, data$theta, label="sgd", type="mse-param") +
  geom_hline(yintercept=1.5, color="green")