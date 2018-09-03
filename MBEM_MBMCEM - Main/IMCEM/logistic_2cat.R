# setwd("/Users/karimimohammedbelhal/Desktop/package_contrib/saemixB/R")
# load("2cat_noavg.RData")
# load("final_cat.RData")
# save.image("final_cat.RData")

setwd("/Users/karimimohammedbelhal/Desktop/imcem_logistic/R")
  source('aaa_generics.R') 
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_distcond.R') 
  source('func_FIM.R')
  source('func_plots.R') 
  source('func_simulations.R') 
  source('main.R')
  source('main_estep.R')
  source('main_estep_incremental.R')
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 

library("mlxR")
library("tensorflow")
library("sgmcmc")
require(ggplot2)
require(gridExtra)
require(reshape2)

setwd("/Users/karimimohammedbelhal/Desktop/imcem_logistic")
source('plots.R') 

# dataset = list("x" = rnorm(1000))
# params = list("theta" = 0)
# logLik = function(params, dataset) {
# distn = tf$contrib$distributions$Normal(params$theta, 1)
# return(tf$reduce_sum(distn$log_prob(dataset$x)))
# }
# stepsize = list("theta" = 1e-4)
# output = sgld(logLik, dataset, params, stepsize)



#####################################################################################
cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Desktop/imcem_logistic/data/logistic_imcem.csv", header=T, sep=",")
# cat_data.saemix<-cat_data.saemix[c(1:1500,6001:7500,12001:13500),]
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"), name.predictors=c("y","dose","time"))


cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]
dose<-xidep[,2]
time<-xidep[,3]

beta0 <- psi[id,1]
gamma0 <- psi[id,2]
delta0 <- psi[id,3]

lm0 <- beta0+gamma0*time + delta0*dose

D <- exp(lm0)+1
P0 <- exp(lm0)/D
P1 <- 1/D

P.obs = (level==0)*P0+(level==1)*P1

return(P.obs)
}

cov <- matrix(c(1,0,0,
                0,1,0,
                0,0,1),ncol=3, byrow=TRUE)

saemix.model<-saemixModel(model=cat_data.model,description="cat model",  type="likelihood" ,
  psi0=matrix(c(2,1,2),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("beta0",
    "gamma0",
    "delta0"))), 
  transform.par=c(0,0,0), fixed.estim=c(1,1,1),covariance.model=cov,omega.init=cov)


K1 =200
K2 = 1

iterations = 1:(K1+K2+1)
end = K1+K2
seed0 = 39546





options<-list(seed=seed0,map=F,fim=F,save.graphs=FALSE, nbiter.mcmc = c(0,0,0,0,1),ll.is=F,displayProgress=TRUE,nb.chains = 1,
              nbiter.saemix = c(K1,K2), map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=100,sampling="seq")
theo.polya<-data.frame(saemix(saemix.model,saemix.data,options))
theo.polya <- cbind(iterations, theo.polya)
theo.polya[end,]


options<-list(seed=seed0,map=F,fim=F,save.graphs=FALSE, nbiter.mcmc = c(2,2,200,0,0),ll.is=F,displayProgress=TRUE,nb.chains = 1,
              nbiter.saemix = c(K1,K2), map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=100,sampling="seq")
theo_ref2<-data.frame(saemix(saemix.model,saemix.data,options))
theo_ref2 <- cbind(iterations, theo_ref2)
theo_ref2[end,]


options.incremental<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(10,20,0,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=50,sampling='randompass')
theo_mix<-data.frame(saemix(saemix.model,saemix.data,options.incremental))
theo_mix <- cbind(iterations, theo_mix)
theo_mix[end,]


options.incremental75<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, nbiter.mcmc = c(10,20,0,0), 
                          nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=75,sampling='randompass')
theo_mix75<-data.frame(saemix(saemix.model,saemix.data,options.incremental75))
theo_mix75 <- cbind(iterations, theo_mix75)


options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(10,100,0,0), nbiter.saemix = c(3*K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=25,sampling='randompass')
theo_mix25<-data.frame(saemix(saemix.model,saemix.data,options.incremental25))
theo_mix25 <- cbind(iterations, theo_mix25)
theo_mix25 <- theo_mix25bis

# theo_mix25 <- theo_mix25second
graphConvMC_3(theo_ref[1:K1,],theo_ref[1:K1,],theo_mix25[1:K1,])
graphConvMC_3(theo_ref[1:K1,],theo_ref[1:K1,],theo_mix25_scaled[1:K1,])
graphConvMC_3(theo_mix25_scaled,theo_mix25_scaled,theo_mix25_scaled)
graphConvMC_3(theo_mix25,theo_mix25,theo_mix25)
graphConvMC_3(theo_ref,theo_ref,theo_ref)

options.incremental40<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(10,20,0,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=40,sampling='randompass')
theo_mix40<-data.frame(saemix(saemix.model,saemix.data,options.incremental40))
theo_mix40 <- cbind(iterations, theo_mix40)


graphConvMC_5(theo_ref,theo_mix25bis,theo_mix40,theo_mix,theo_mix75)
graphConvMC_3(theo_ref,theo_mix75,theo_mix40)

theo_ref_scaled <- theo_ref
theo_mix50_scaled <- theo_mix
theo_mix75_scaled <- theo_mix75
theo_mix25_scaled <- theo_mix25
theo_mix40_scaled <- theo_mix40


theo_ref_scaled$iterations = theo_ref_scaled$iterations*1
theo_mix50_scaled$iterations = theo_mix50_scaled$iterations*0.5
theo_mix75_scaled$iterations = theo_mix75_scaled$iterations*0.75
theo_mix25_scaled$iterations = theo_mix25_scaled$iterations*0.25
theo_mix40_scaled$iterations = theo_mix40_scaled$iterations*0.4
graphConvMC_3(theo_ref_scaled,theo_mix50_scaled,theo_mix25_scaled)

graphConvMC_3(theo_ref_scaled,theo_mix50_scaled,theo_mix40_scaled)
graphConvMC_5(theo_ref_scaled,theo_mix25_scaled,theo_mix40_scaled,theo_mix50_scaled,theo_mix75_scaled)


options.incremental85<-list(seed=seed0,map=F,fim=F,ll.is=F,save.graphs=FALSE,nb.chains = 1, 
  nbiter.mcmc = c(10,20,0,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),
  nbiter.sa=0,nbiter.burn =0, nb.replacement=85)
theo_mix85<-data.frame(saemix(saemix.model,saemix.data,options.incremental85))
theo_mix85 <- cbind(iterations, theo_mix85)


theo_mix85_scaled <- theo_mix85
theo_mix85_scaled$iterations = theo_mix85_scaled$iterations*0.85
graphConvMC_6(theo_ref_scaled[10:end,],
  theo_mix25_scaled[10:end,],theo_mix40_scaled[10:end,],
  theo_mix50_scaled[10:end,],theo_mix75_scaled[10:end,],
  theo_mix85_scaled[10:end,])



seplot <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  graf <- ggplot(df)+geom_line(aes(iterations,value,by=value,colour = df$algo),show.legend = legend) +
  xlab("iterations") + ylab('value') + facet_wrap(~variable,scales = "free_y") + theme_bw() 
  grid.arrange(graf)
  # do.call("grid.arrange", c(graf, ncol=1, top=title))
}

theo_ref_scaledbis <- theo_ref_scaled
theo_mix50_scaledbis <- theo_mix50_scaled
theo_mix75_scaledbis <- theo_mix75_scaled
theo_mix25_scaledbis <- theo_mix25_scaled
theo_mix40_scaledbis <- theo_mix40_scaled
theo_mix85_scaledbis <- theo_mix85_scaled

theo_ref_scaledbis$algo = 'ref'
theo_mix50_scaledbis$algo = '50'
theo_mix75_scaledbis$algo = '75'
theo_mix25_scaledbis$algo = '25'
theo_mix40_scaledbis$algo = '40'
theo_mix85_scaledbis$algo = '85'



for (i in 2:4){
comparison <- 0
comparison <- rbind(theo_ref_scaledbis[10:end,c(1,i,8)],
  theo_mix25_scaledbis[10:end,c(1,i,8)],theo_mix40_scaledbis[10:end,c(1,i,8)],
  theo_mix50_scaledbis[10:end,c(1,i,8)],theo_mix75_scaledbis[10:end,c(1,i,8)],
  theo_mix85_scaledbis[10:end,c(1,i,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)


prec <- seplot(var, title="comparison",legend=TRUE)
assign(paste("prec", i, sep = ""), prec) 
# setwd("/Users/karimimohammedbelhal/Desktop/")
# ggsave(paste("precwarfa_seq_50sim_100indiv_", i, ".png", sep=""),prec)
}

grid.arrange(prec2,prec3,prec4, ncol=2)

