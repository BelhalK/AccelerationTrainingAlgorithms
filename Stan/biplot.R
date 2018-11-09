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
# load("RData/hmc_quantile_indiv.RData")
library(mvtnorm)
library(MASS)
# save.image("biplotviz.RData")
load("biplotviz.RData")
source("marginal_plot.R")
#### PLOT OF THE TRUTH (CURVE) #####

# yobs = warfa_data[warfa_data$id==11,4]
# xreg = warfa_data[warfa_data$id==11,2:3]
# model1cpt<-function(psi,xidep) { 
#   time<-xidep[,1]
#   dose<-xidep[,2]
#   ka<-psi[1]
#   V<-psi[2]
#   k<-psi[3]

#   ypred<-dose*ka/(V*(ka-k))*(exp(-k*time)-exp(-ka*time))
#   return(ypred)
# }

# fpred <- model1cpt(etamap[10,],xreg)

# sigma.prior <- matrix(c(sqrt(0.2),0,0,0,sqrt(0.18),0,0,0,sqrt(0.03)),ncol=3,byrow=TRUE)
# sigma.error <- diag(11)

# posterior <- function(psi) { 
#   dens <- dmvnorm(yobs-model1cpt(psi,xreg), mean = rep(0,11),sigma = sigma.error)*dmvnorm(psi, mean =  rep(0,3), sigma = sigma.prior)
#   return(dens)
# }

# exp(-norm(as.matrix(yobs-model1cpt(etamap[1,],xreg))))
# dmvnorm(etamap[1,], mean =  rep(0,3), sigma = sigma.prior)

# posterior(etamap[1,])



samples.imh.proposal <- data.frame(mvrnorm(n = 10000000, map, G.map, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
samples.advi.proposal <- data.frame(mvrnorm(n = 10000000, mu.vi, Gamma.vi, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
head(vi)

samples.imh.proposal$algo <- "IMH"
samples.advi.proposal$algo <- "VI"
vi$algo <- "Truth"

colnames(samples.imh.proposal) <- colnames(samples.advi.proposal) <- colnames(vi) <- c("ka","V","k","algo")



d.newka <- density(samples.imh.proposal$ka) # returns the density data 
d.advika <- density(samples.advi.proposal$ka) # returns the density data 
d.vika <- density(vi$ka) # returns the density data 


d.newV <- density(samples.imh.proposal$V) # returns the density data 
d.adviV <- density(samples.advi.proposal$V) # returns the density data 
d.viV <- density(vi$V) # returns the density data 


d.newk <- density(samples.imh.proposal$k) # returns the density data 
d.advik <- density(samples.advi.proposal$k) # returns the density data 
d.vik <- density(vi$k) # returns the density data 


par(mfrow=c(1,3))
plot(d.advika,col="blue", main = "ka", xlab="")
lines(d.newka,col="red")
lines(d.vika,col="green")
# legend('topright', c("VI","IMH","Posterior") , lty=1, col=c("blue","red","yellow"), bty='n', cex=1)

plot(d.adviV,col="blue", main = "V", xlab="",ylab="")
lines(d.newV,col="red")
lines(d.viV,col="green")
# legend('topright', c("VI","IMH","Posterior") , lty=1, col=c("blue","red","green"), bty='n', cex=1)

plot(d.advik,col="blue", main = "k", xlab="",ylab="")
lines(d.newk,col="red")
lines(d.vik,col="green")
legend('topright', c("VI","IMH","Posterior") , lty=1, col=c("blue","red","green"), bty='n', cex=1)


library(ggstatsplot)

ggscatterstats(
  data = iris,                                          
  x = Sepal.Length,                                                  
  y = Sepal.Width,
  xlab = "Sepal Length",
  ylab = "Sepal Width",
  marginal = TRUE,
  marginal.type = "histogram",
  centrality.para = "mean",
  margins = "both",
  title = "Relationship between Sepal Length and Sepal Width",
  messages = FALSE
)


data <- rbind(samples.imh.proposal[1:1000000,],samples.advi.proposal[1:1000000,],vi[1:20000,])
# data <- rbind(samples.imh.proposal[1:100000,],samples.advi.proposal[1:100000,],vi[1:20000,])

colnames(data) <- c("ka","V","k","Proposal")


marginal_plot(x = ka, y = V, 
  group = Proposal, data = data, bw = "nrd", 
  lm_formula = NULL, xlab = "ka", ylab = "V", pch = 15, cex = 0.5)
dev.copy(jpeg,'biplotkaV.jpg', width=900, height=550)
dev.off()


marginal_plot(x = ka, y = k, 
  group = Proposal, data = data, bw = "nrd", 
  lm_formula = NULL, xlab = "ka", ylab = "k", pch = 15, cex = 0.5)
dev.copy(jpeg,'biplotkak.jpg', width=900, height=550)
dev.off()

marginal_plot(x = V, y = k, 
  group = Proposal, data = data, bw = "nrd", 
  lm_formula = NULL, xlab = "V", ylab = "k", pch = 15, cex = 0.5)
dev.copy(jpeg,'biplotkV.jpg', width=900, height=550)
dev.off()





#### PLOTS OF THE PROPOSALS AGAINST THE TRUTH #####

# dmvnorm(x=c(0,0))
# dmvnorm(x=c(0,0), mean=c(1,1))
# x <- rmvnorm(n=100, mean=c(1,1))
# plot(x)

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
lines(d.vi) # plo

# d.vi <- density(vi$ka/0.002) # returns the density data 
# plot(d.vi) # plo


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




sigma.imh <- matrix(c(sqrt(G.map[2,2]),0,0,sqrt(G.map[3,3])), nrow=2)
sigma.advi <- matrix(c(sqrt(Gamma.vi[2,2]),0,0,sqrt(Gamma.vi[3,3])), nrow=2)
data.grid <- expand.grid(s.1 = seq(-1, 1, length.out=200), s.2 = seq(-1, 1, length.out=200))
q.imh <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean =  map[c(2,3)], sigma = sigma.imh))
q.advi <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu.vi[c(2,3)], sigma = sigma.advi))

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