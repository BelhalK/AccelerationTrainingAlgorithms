#library(rstan)
setwd("/Users/karimimohammedbelhal/Desktop/variationalBayes/mcmc_R_isolate/Dir2")
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_cov.R') 
  source('func_distcond.R') 
  source('func_FIM.R') 
  source('func_ggplot2.R') 
  source('func_plots.R') 
  source('func_simulations.R') 
  source('ggplot2_global.R') 
  # source('KL.R') 
  #source('vi.R') 
  source('global.R')
  source('main.R')
  source('mcmc_main.R') 
  source('main_estep.R')
  source('main_estep_mcmc.R') 
  source('main_estep_morekernels.R') 
  # source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('plots_ggplot2.R') 
  source('saemix-package.R') 
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  source("mixtureFunctions.R")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mcmc_newkernel")
source('mcmc.R')
source('mcmc_mix.R')
source('mcmc_sum.R')
source('initalgo.R') 
source("mixtureFunctions.R")
library("mlxR")

require(ggplot2)
require(gridExtra)
require(reshape2)
library(grid)
library(lattice)
#####################################################################################


# Doc
# oxboys.saemix<-read.table( "data/oxboys.saemix.tab",header=T,na=".")
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/oxboys")
oxboys.saemix<-read.table( "ox_synth.csv",header=T,na=".",sep=",")
oxboys.saemix_less <- oxboys.saemix[1:10,1:3]
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/mcmc_newkernel")
saemix.data<-saemixData(name.data=oxboys.saemix_less,header=TRUE,
  name.group=c("id"),name.predictors=c("time"),name.response=c("y"),
  units=list(x="yr",y="cm"))


growth.linear<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (2 columns, base and slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  base<-psi[id,1]
  slope<-psi[id,2]
  f<-base+slope*x
  return(f)
}
saemix.model<-saemixModel(model=growth.linear,description="Linear model",
  psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(1,0),covariance.model=matrix(c(1,1,1,1),ncol=2,byrow=TRUE), 
  error.model="constant")




indiv = 1
seed0 = 35644
replicate = 50
iter_mcmc = 1000
burn = 400


saemix.options_rwm<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,iter_mcmc,iter_mcmc,0))
saemix.options_linear<-list(seed=seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))

#reference rwm
ref <- mcmc(saemix.model,saemix.data,saemix.options_rwm,iter_mcmc)
new<-mcmc(saemix.model,saemix.data,saemix.options_linear,iter_mcmc)

graphConvMC_twokernels(new$eta[[indiv]],ref$eta[[indiv]], title="eta")

final_rwm <- 0
expec_rwm <- ref$eta[[indiv]]
var_rwm <- ref$eta[[indiv]]
expec_rwm[,2:3] <- 0 
var_rwm[,2:3] <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_rwm<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,iter_mcmc,iter_mcmc,0))
  post_rwm<-mcmc(saemix.model,saemix.data,saemix.options_rwm)$eta[[indiv]]
  post_rwm['individual'] <- j
  expec_rwm[,2:3] <- expec_rwm[,2:3] + post_rwm[,2:3]
  var_rwm[,2] <- var_rwm[,2] + (post_rwm[,2])^2
  var_rwm[,3] <- var_rwm[,3] + (post_rwm[,3])^2
  final_rwm <- rbind(final_rwm,post_rwm)
}
expec_rwm[,2:3] <- expec_rwm[,2:3]/replicate
var_rwm[,2:3] <- var_rwm[,2:3]/replicate



names(final_rwm)[1]<-paste("time")
names(final_rwm)[4]<-paste("id")
final_rwm1 <- final_rwm[c(4,1,2)]
final_rwm2 <- final_rwm[c(4,1,3)]



# prctilemlx(final_rwm1[-1,],band = list(number = 8, level = 80))
# prctilemlx(final_rwm2[-1,],band = list(number = 8, level = 80))
# prctilemlx(final_rwm3[-1,],band = list(number = 8, level = 80))

final_new <- 0
expec_new <- new$eta[[indiv]]
var_new <- new$eta[[indiv]]
expec_new[,2:3] <- 0 
var_new[,2:3] <- 0
for (j in 1:replicate){
  print(j)
  saemix.options_linear<-list(seed=j*seed0,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(0,0,0,iter_mcmc))
  post_new<-mcmc(saemix.model,saemix.data,saemix.options_linear)$eta[[indiv]]
  post_new['individual'] <- j
  expec_new[,2:3] <- expec_new[,2:3] + post_new[,2:3]
  var_new[,2] <- var_new[,2] + (post_new[,2])^2
  var_new[,3] <- var_new[,3] + (post_new[,3])^2
  final_new <- rbind(final_new,post_new)
}
expec_new[,2:3] <- expec_new[,2:3]/replicate
var_new[,2:3] <- var_new[,2:3]/replicate




names(final_new)[1]<-paste("time")
names(final_new)[4]<-paste("id")
final_new1 <- final_new[c(4,1,2)]
final_new2 <- final_new[c(4,1,3)]



final_rwm1['group'] <- 1
final_new1['group'] <- 2
final_new1$id <- final_new1$id +1

final1 <- 0
final1 <- rbind(final_rwm1[-1,],final_new1[-1,])
labels <- c("ref","new")
# prctilemlx(final1[c(1,4,2,3)], band = list(number = 4, level = 80),group='group', label = labels) 
# plt1 <- prctilemlx(final1, band = list(number = 4, level = 80),group='group', label = labels) 

# rownames(final1) <- 1:nrow(final1)

plot.S1 <- plot.prediction.intervals(final1[c(1,4,2,3)], 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('#01b7a5', '#c17b01'))
plot.S <- plot.S1  + ylab("base")+ theme(legend.position=c(0.9,0.8))+ theme_bw()
# print(plot.S1)



final_rwm2['group'] <- 1
final_new2['group'] <- 2
final_new2$id <- final_new2$id +1

final2 <- 0
final2 <- rbind(final_rwm2[-1,],final_new2[-1,])
labels <- c("ref","new")
# prctilemlx(final2[c(1,4,2,3)], band = list(number = 4, level = 80),group='group', label = labels) 
# plt1 <- prctilemlx(final1, band = list(number = 4, level = 80),group='group', label = labels) 

# rownames(final1) <- 1:nrow(final1)

plot.S2 <- plot.prediction.intervals(final2[c(1,4,2,3)], 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('#01b7a5', '#c17b01'))
plot.S2 <- plot.S2  + ylab("slope")+ theme(legend.position=c(0.9,0.8))+ theme_bw()

grid.arrange(plot.S, plot.S2,ncol=2)


graphConvMC_twokernels(expec_rwm,expec_new, title="Expectations")
graphConvMC_twokernels(var_rwm,var_new, title="Variances")




#values table
post_mean_rwm <- 0
var_rwm <- 0
error_rwm <- 0
true_param <- colMeans(final_rwm[burn:iter_mcmc,2:3])
for (j in 1:replicate){
  post_mean_rwm <- post_mean_rwm + colMeans(final_rwm[(j*burn):(j*iter_mcmc),2:3])
}
post_mean_rwm = 1/replicate*post_mean_rwm

for (j in 1:replicate){
  var_rwm <- var_rwm + (final_rwm[(j*iter_mcmc),2:3]-post_mean_rwm)^2
  error_rwm <- error_rwm + (final_rwm[(j*iter_mcmc),2:3]-true_param)^2
}

error_rwm = 1/replicate*error_rwm
var_rwm = 1/replicate*var_rwm



#values table
post_mean_mix <- 0
var_mix <- 0
error_mix <- 0
true_param <- colMeans(final_rwm[burn:iter_mcmc,2:3])
for (j in 1:replicate){
  post_mean_mix <- post_mean_mix + colMeans(final_new[(j*burn):(j*iter_mcmc),2:3])
}
post_mean_mix = 1/replicate*post_mean_mix

for (j in 1:replicate){
  var_mix <- var_mix + (final_new[(j*iter_mcmc),2:3]-post_mean_mix)^2
  error_mix <- error_mix + (final_new[(j*iter_mcmc),2:3]-true_param)^2
}

error_mix = 1/replicate*error_mix
var_mix = 1/replicate*var_mix

