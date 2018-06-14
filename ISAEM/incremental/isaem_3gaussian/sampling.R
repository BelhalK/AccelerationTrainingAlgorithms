
# invcdf <- function(u, m) {
#     return(1+exp(-u*m))
# }

invcdf <- function(u, m) {
    return(sqrt(m^2/(1 - (1 - m^2) * u)))
}


sample1 <- sapply(runif(100), invcdf, m = 0.2)
plot(density(sample1), main = "Sample Density using invcdf Function")



for (i in (1:length(sample1)))
{
	if (sample1[i] < 1)
		print(sample1[i])
}

l <- rmultinom(n=10,size=99,prob=c(sample1[1],1-sample1[1]))
length(l[1,])


endsign <- function(f, sign = 1) {
    b <- sign
    while (sign * f(b) < 0) b <- 10 * b
    return(b)
}

samplepdf <- function(n, pdf, ..., spdf.lower = -Inf, spdf.upper = Inf) {
    vpdf <- function(v) sapply(v, pdf, ...)  # vectorize
    cdf <- function(x) integrate(vpdf, spdf.lower, x)$value
    invcdf <- function(u) {
        subcdf <- function(t) cdf(t) - u
        if (spdf.lower == -Inf) 
            spdf.lower <- endsign(subcdf, -1)
        if (spdf.upper == Inf) 
            spdf.upper <- endsign(subcdf)
        return(uniroot(subcdf, c(spdf.lower, spdf.upper))$root)
    }
    sapply(runif(n), invcdf)
}

h <- function(t, m) {
    if (t >= m & t <= 1) 
        return(2 * m^2/(1 - m^2)/t^3)
    return(0)
}

sigmoid <- function(t, m) {
    if (t >= m & t <= 1) 
        return(1/(1+exp(-t)))
    return(0)
}

sample2 <- samplepdf(100, h, m = .5)
a <- samplepdf(10, h, m = 0.3)
samplepdf(1, sigmoid, m = 0.1)


rmultinom(n=1,size=99,prob=c(a[1],1-a[1]))