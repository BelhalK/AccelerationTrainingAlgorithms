T=1e5
prop=x=rnorm(T,sd=.01)
ratop=dnorm(prop,log=TRUE)-dnorm(prop,sd=.01,log=TRUE)
ratav=ratop[1]
logu=ratop-log(runif(T))
for (t in 2:T){
  if (logu[t]>ratav){
    x[t]=prop[t];ratav=ratop[t]}else{x[t]=x[t-1]}
  }

plot(x)
