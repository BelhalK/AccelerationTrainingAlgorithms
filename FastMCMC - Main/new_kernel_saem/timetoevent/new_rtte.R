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
  # source('SaemixObject.R') 
  source('zzz.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent")
source('main_time.R')
source('main_estep_time.R')
source('main_mstep_time.R') 
source('func_aux_time.R') 
source('SaemixObject_time.R') 
source('main_initialiseMainAlgo_time.R') 
source("mixtureFunctions.R")

library("mlxR")
library("psych")
library("coda")
library("Matrix")
library(abind)
require(ggplot2)
require(gridExtra)
require(reshape2)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")


timetoevent.saemix <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent/rtte1.csv", header=T, sep=",")
timetoevent.saemix <- timetoevent.saemix[timetoevent.saemix$ytype==2,]
# timetoevent.saemix["nb"] <- 0
# for (i in 1:length(unique(timetoevent.saemix$id))) {
#     timetoevent.saemix[timetoevent.saemix$id==i,5] <- length(which(timetoevent.saemix[timetoevent.saemix$id==i,3]==1))
#   }

saemix.data<-saemixData(name.data=timetoevent.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("time","y"), name.X=c("time"))
# write.table(timetoevent.saemix[,1:3],"rtte.txt",sep=",",row.names=FALSE)


timetoevent.model<-function(psi,id,xidep) {
T<-xidep[,1]
y<-xidep[,2]
N <- nrow(psi)
Nj <- length(T)

censoringtime = 20

lambda <- psi[id,1]
beta <- psi[id,2]

init <- which(T==0)
cens <- which(T==censoringtime)
ind <- setdiff(1:Nj, append(init,cens))


hazard <- (beta/lambda)*(T/lambda)^(beta-1)
H <- (T/lambda)^beta

logpdf <- rep(0,Nj)
logpdf[cens] <- -H[cens] + H[cens-1]
logpdf[ind] <- -H[ind] + H[ind-1] + log(hazard[ind])

return(logpdf)
}


saemix.model<-saemixModel(model=timetoevent.model,description="time model",   
  psi0=matrix(c(10,1),ncol=2,byrow=TRUE,dimnames=list(NULL,   
  c("lambda","beta"))), 
  transform.par=c(1,1),covariance.model=matrix(c(1,0,0,1),ncol=2, 
  byrow=TRUE))


K1 = 450
K2 = 50

iterations = 1:(K1+K2+1)
gd_step = 0.01
end = K1+K2
#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE, map.range=c(0))
theo_ref<-data.frame(saemix_time(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

# graphConvMC_saem(theo_ref, title="new kernel")

#ref (map always)
options.cat<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(1:10))
cat_saem<-data.frame(saemix_time(saemix.model,saemix.data,options.cat))
cat_saem <- cbind(iterations, cat_saem)

# graphConvMC_saem(cat_saem, title="new kernel")
graphConvMC2_saem(theo_ref,cat_saem, title="new kernel")

seed0 = 39546
replicate = 10
final_rwm <- 0
for (j in 1:replicate){
  print(j)
  options<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=TRUE, map.range=c(0))
  theo_ref<-data.frame(saemix_time(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- j
  final_rwm <- rbind(final_rwm,theo_ref)
}

names(final_rwm)[1]<-paste("time")
names(final_rwm)[6]<-paste("id")
final_rwm1 <- final_rwm[c(6,1,2)]
final_rwm2 <- final_rwm[c(6,1,3)]
final_rwm3 <- final_rwm[c(6,1,4)]
final_rwm4 <- final_rwm[c(6,1,5)]

# prctilemlx(final_rwm1[-1,],band = list(number = 8, level = 80)) + ggtitle("RWM")

#mix (RWM and MAP new kernel for liste of saem iterations)
final_mix <- 0
for (j in 1:replicate){
  print(j)
  options.mix<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,6),nbiter.saemix = c(K1,K2),displayProgress=FALSE, map.range=c(1:10))
  theo_mix<-data.frame(saemix_time(saemix.model,saemix.data,options.mix))
  theo_mix <- cbind(iterations, theo_mix)
  theo_mix['individual'] <- j
  final_mix <- rbind(final_mix,theo_mix)
}


names(final_mix)[1]<-paste("time")
names(final_mix)[6]<-paste("id")
final_mix1 <- final_mix[c(6,1,2)]
final_mix2 <- final_mix[c(6,1,3)]
final_mix3 <- final_mix[c(6,1,4)]
final_mix4 <- final_mix[c(6,1,5)]


# prctilemlx(final_mix1[-1,1:3],band = list(number = 8, level = 80)) + ggtitle("mix")






final_rwm1['group'] <- 1
final_mix1['group'] <- 2
final_mix1$id <- final_mix1$id +1


final1 <- rbind(final_rwm1[-1,],final_mix1[-1,])
labels <- c("ref","new")
# prctilemlx(final1[c(1,4,2,3)], band = list(number = 4, level = 80),group='group', label = labels) 
# plt1 <- prctilemlx(final1, band = list(number = 4, level = 80),group='group', label = labels) 

# rownames(final1) <- 1:nrow(final1)

plot.S1 <- plot.prediction.intervals(final1[c(1,4,2,3)], 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('#01b7a5', '#c17b01'))
plot.S <- plot.S1  + ylab("lambda")+ theme(legend.position=c(0.9,0.8))+ theme_bw()
# print(plot.S1)



final_rwm2['group'] <- 1
final_mix2['group'] <- 2
final_mix2$id <- final_mix2$id +1


final2 <- rbind(final_rwm2[-1,],final_mix2[-1,])
labels <- c("ref","new")
# prctilemlx(final2[c(1,4,2,3)], band = list(number = 4, level = 80),group='group', label = labels) 
# plt1 <- prctilemlx(final1, band = list(number = 4, level = 80),group='group', label = labels) 

# rownames(final1) <- 1:nrow(final1)

plot.S2 <- plot.prediction.intervals(final2[c(1,4,2,3)], 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('#01b7a5', '#c17b01'))
plot.S2 <- plot.S2  + ylab("beta")+ theme(legend.position=c(0.9,0.8))+ theme_bw()


final_rwm3['group'] <- 1
final_mix3['group'] <- 2
final_mix3$id <- final_mix3$id +1


final3 <- rbind(final_rwm3[-1,],final_mix3[-1,])
labels <- c("ref","new")
# prctilemlx(final3[c(1,4,2,3)], band = list(number = 4, level = 80),group='group', label = labels) 
# plt1 <- prctilemlx(final1, band = list(number = 4, level = 80),group='group', label = labels) 

# rownames(final1) <- 1:nrow(final1)

plot.S3 <- plot.prediction.intervals(final3[c(1,4,2,3)], 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('#01b7a5', '#c17b01'))
plot.S3 <- plot.S3  + ylab("w2_lambda")+ theme(legend.position=c(0.9,0.8))+ theme_bw()




final_rwm4['group'] <- 1
final_mix4['group'] <- 2
final_mix4$id <- final_mix4$id +1


final4 <- rbind(final_rwm4[-1,],final_mix4[-1,])
labels <- c("ref","new")
# prctilemlx(final4[c(1,4,2,3)], band = list(number = 4, level = 80),group='group', label = labels) 
# plt1 <- prctilemlx(final1, band = list(number = 4, level = 80),group='group', label = labels) 

# rownames(final1) <- 1:nrow(final1)

plot.S4 <- plot.prediction.intervals(final4[c(1,4,2,3)], 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('#01b7a5', '#c17b01'))
plot.S4 <- plot.S4  + ylab("w2_beta")+ theme(legend.position=c(0.9,0.8))+ theme_bw()



grid.arrange(plot.S, plot.S2,plot.S3,plot.S4,ncol=3)


#values table

#values table
sample_mean_rwm <- 0
var_rwm <- 0
error_rwm <- 0
ka_true =  1.28289
V_true = 7.9396
k_true = 0.13452/7.9396

o_ka_true =  0.73671^2
o_V_true = 0.13672^2
o_k_true = 0.26209^2
a_true = 0.26629
replicate = 10

true_param <- c(ka_true,V_true,k_true,o_ka_true,o_V_true,o_k_true,a_true)
for (j in 1:replicate){
  sample_mean_rwm <- sample_mean_rwm + colMeans(final_rwm[(j*K1):(j*(K1+K2)),c(2,3,4,5,6,7,8)])
}
sample_mean_rwm = 1/replicate*sample_mean_rwm

for (j in 1:replicate){
  var_rwm <- var_rwm + (final_rwm[(j*(K1+K2)),c(2,3,4,5,6,7,8)]-sample_mean_rwm)^2
  error_rwm <- error_rwm + (final_rwm[(j*(K1+K2)),c(2,3,4,5,6,7,8)]-true_param)^2
}

error_rwm = 1/replicate*error_rwm
var_rwm = 1/replicate*var_rwm




sample_mean_mix <- 0
var_mix <- 0
error_mix <- 0

for (j in 1:replicate){
  sample_mean_mix <- sample_mean_mix + colMeans(final_mix[(j*K1):(j*(K1+K2)),c(2,3,4,5,6,7,8)])
}
sample_mean_mix = 1/replicate*sample_mean_mix

for (j in 1:replicate){
  var_mix <- var_mix + (final_mix[(j*(K1+K2)),c(2,3,4,5,6,7,8)]-sample_mean_mix)^2
  error_mix <- error_mix + (final_mix[(j*(K1+K2)),c(2,3,4,5,6,7,8)]-true_param)^2
}

error_mix = 1/replicate*error_mix
var_mix = 1/replicate*var_mix






plot.prediction.intervals <- function(r, plot.median=TRUE, level=1, labels=NULL, 
                                      legend.title=NULL, colors=NULL) {
  P <- prctilemlx(r, number=1, level=level, plot=FALSE)
  if (is.null(labels))  labels <- levels(r$group)
  if (is.null(legend.title))  legend.title <- "group"
  names(P$y)[2:4] <- c("p.min","p50","p.max")
  pp <- ggplot(data=P$y)+ylab(NULL)+ 
    geom_ribbon(aes(x=time,ymin=p.min, ymax=p.max,fill=group),alpha=.5) 
  if (plot.median)
    pp <- pp + geom_line(aes(x=time,y=p50,colour=group))
  
  if (is.null(colors)) {
    pp <- pp + scale_fill_discrete(name=legend.title,
                                   breaks=levels(r$group),
                                   labels=labels)
    pp <- pp + scale_colour_discrete(name=legend.title,
                                     breaks=levels(r$group),
                                     labels=labels, 
                                     guide=FALSE)
  } else {
    pp <- pp + scale_fill_manual(name=legend.title,
                                 breaks=levels(r$group),
                                 labels=labels,
                                 values=colors)
    pp <- pp + scale_colour_manual(name=legend.title,
                                   breaks=levels(r$group),
                                   labels=labels,
                                   guide=FALSE,values=colors)
  }  
  return(pp)
}

