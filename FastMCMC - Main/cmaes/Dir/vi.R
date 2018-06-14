vi <- function() {
    stanmodelcode <- "
	data {
	  int<lower=0> N;
	  int<lower=0> y[N];
	}
	parameters {
	  real<lower=0.00001> Theta;
	}
	model {
	  Theta ~ normal(0, 1);
	  for (n in 1:N)
	    y[n] ~ exponential(Theta);
	}
	generated quantities{
	  real y_pred;
	  y_pred = exponential_rng(Theta);
	}
	"
	stanDso <- stan_model(model_code = stanmodelcode) 
	claims.obs <- c(100, 950, 450)
	N <- length(claims.obs)
	dat <- list(N = N, y = claims.obs); 

	fit <- sampling(stanDso, data = dat, iter = 10000, warmup=200) 

	Theta <- extract(fit, 'Theta')
	Theta <- unlist(Theta, use.names=FALSE)
	y_pred <- extract(fit, 'y_pred')
	y_pred <- unlist(y_pred, use.names=FALSE)
	return(list(ypred = y_pred, theta = Theta))
    }
