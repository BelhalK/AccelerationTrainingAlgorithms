library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)
library(dplyr)
library(data.table)
library(rstan)
load("RData/hmc_quantile_indiv.RData")
library(mvtnorm)

#### PLOT OF THE TRUTH (CURVE) #####

yobs = warfa_data[warfa_data$id==11,4]
x = warfa_data[warfa_data$id==11,2:3]
model1cpt<-function(psi,xidep) { 
  time<-xidep[,1]
  dose<-xidep[,2]
  ka<-psi[1]
  V<-psi[2]
  k<-psi[3]

  ypred<-dose*ka/(V*(ka-k))*(exp(-k*time)-exp(-ka*time))
  return(ypred)
}

fpred <- model1cpt(etamap[10,],x)

sigma.prior <- matrix(c(sqrt(0.2),0,0,0,sqrt(0.18),0,0,0,sqrt(0.03)),ncol=3,byrow=TRUE)
sigma.error <- daig(11)
posterior <- function(psi) { 
  dens <- dmvnorm(yobs-model1cpt(psi,x), rep(0,11),sigma = sigma.error)*dmvnorm(psi, mean =  rep(0,3), sigma = sigma.prior)
  
  return(dens)
}


posterior(etamap[1,])




#### PLOTS OF THE PROPOSALS AGAINST THE TRUTH #####

dmvnorm(x=c(0,0))
dmvnorm(x=c(0,0), mean=c(1,1))
x <- rmvnorm(n=100, mean=c(1,1))
plot(x)

mu.vi 
map <- etamap[10,] 

Gamma.vi 
G.map <- Gammamap[[10]]

colnames(mu.vi) <- colnames(map) <- colnames(Gamma.vi) <- colnames(G.map) <- c("ka","V","k")

plot(d.new)
lines(d.advi)
lines(d.vi)




p1 <- ggplot(data = data.frame(x = c(-1, 1)), aes(x)) +
  stat_function(fun = dnorm, n = 500, args = list(mean = mu.vi[1], sd = sqrt(Gamma.vi[1,1])), colour="red")+
  stat_function(fun = dnorm, n = 500, args = list(mean = map[1], sd = sqrt(G.map[1,1])), colour="blue") + 
  stat_function(fun=function(x) 1+ 3 * x, colour = "yellow") +
  ylab("") +
  scale_y_continuous(breaks = NULL)
p1


# d.vi <- density(vi$ka) # returns the density data 
# lines(d.vi)


p2 <- ggplot(data = data.frame(x = c(-1, 1)), aes(x)) +
  stat_function(fun = dnorm, n = 500, args = list(mean = mu.vi[2], sd = sqrt(Gamma.vi[2,2])), colour="red")+
  stat_function(fun = dnorm, n = 500, args = list(mean = map[2], sd = sqrt(G.map[2,2])), colour="blue") + 
  stat_function(fun=function(x) 1+ 3 * x, colour = "yellow") +
  ylab("") +
  scale_y_continuous(breaks = NULL)
p2
# d.vi <- density(vi$V) # returns the density data 
# lines(d.vi)



p3 <- ggplot(data = data.frame(x = c(-1, 1)), aes(x)) +
  stat_function(fun = dnorm, n = 500, args = list(mean = mu.vi[3], sd = sqrt(Gamma.vi[3,3])), colour="red")+
  stat_function(fun = dnorm, n = 500, args = list(mean = map[3], sd = sqrt(G.map[3,3])), colour="blue") + 
  stat_function(fun=function(x) 1+ 3 * x, colour = "yellow") +
  ylab("") +
  scale_y_continuous(breaks = NULL)
p3
# d.vi <- density(vi$k) # returns the density data 
# lines(d.vi)


sigma.imh <- matrix(c(sqrt(G.map[1,1]),0,0,sqrt(G.map[2,2])), nrow=2)
sigma.advi <- matrix(c(sqrt(Gamma.vi[1,1]),0,0,sqrt(Gamma.vi[2,2])), nrow=2)
data.grid <- expand.grid(s.1 = seq(-1, 1, length.out=200), s.2 = seq(-1, 1, length.out=200))
q.imh <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean =  map[1:2], sigma = sigma.imh))
q.advi <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu.vi[1:2], sigma = sigma.advi))

qq <- rbind(q.imh,q.advi)
qq$algo <- c(rep("IMH",40000),rep("ADVI",40000))

ggplot(qq, aes(x=s.1, y=s.2, z=prob, group = algo, col=algo)) + 
    geom_contour() +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) 



sigma.imh <- matrix(c(sqrt(G.map[1,1]),0,0,sqrt(G.map[3,3])), nrow=2)
sigma.advi <- matrix(c(sqrt(Gamma.vi[1,1]),0,0,sqrt(Gamma.vi[3,3])), nrow=2)
data.grid <- expand.grid(s.1 = seq(-1, 1, length.out=200), s.2 = seq(-1, 1, length.out=200))
q.imh <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean =  map[c(1,3)], sigma = sigma.imh))
q.advi <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu.vi[c(1,3)], sigma = sigma.advi))

qq <- rbind(q.imh,q.advi)
qq$algo <- c(rep("IMH",40000),rep("ADVI",40000))

ggplot(qq, aes(x=s.1, y=s.2, z=prob, group = algo, col=algo)) + 
    geom_contour() +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), ratio = 1) 



########## PLOTS OF THE RESULTING MCMC ######
# imh: new
# variational mcmc: advi
# truth: vi
new$algo <- "IMH"
advi$algo <- "VI"
vi$algo <- "Truth"

colnames(new) <- colnames(advi) <- colnames(vi) <- c("ka","V","k","algo")
d.new <- density(new$ka) # returns the density data 
plot(d.new) # plots the results

d.advi <- density(advi$ka) # returns the density data 
plot(d.advi) # plots the results

d.vi <- density(vi$ka) # returns the density data 
plot(d.vi) # plots the results

plot(d.new)
lines(d.advi)
lines(d.vi)
# Compare MPG distributions for cars with 
# 4,6, or 8 cylinders
library(sm)

# create value labels 
algo.f <- factor(algo, levels= c(4,6,8),
  labels = c("IMH", "VI", "Truth")) 

# plot densities 
sm.density.compare(mpg, cyl, xlab="Miles Per Gallon")
title(main="MPG Distribution by Car Cylinders")

# add legend via mouse click
colfill<-c(2:(2+length(levels(cyl.f)))) 
legend(locator(1), levels(cyl.f), fill=colfill)