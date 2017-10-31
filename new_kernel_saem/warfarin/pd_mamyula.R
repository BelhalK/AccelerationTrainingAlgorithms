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
  source('main_initialiseMainAlgo.R') 
  source('main_mstep.R') 
  source('SaemixData.R')
  source('plots_ggplot2.R') 
  source('saemix-package.R') 
  source('SaemixModel.R') 
  source('SaemixRes.R') 
  source('SaemixObject.R') 
  source('zzz.R') 
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem")
source('newkernel_main.R')
source('main_new.R')
source('main_estep_new.R')
source('main_estep_new2.R')
source('main_gd.R')
source('main_estep_gd.R')
source('main_estep_newkernel.R')
source('main_gd_mix.R')
source('main_estep_gd_mix.R')
source('main_estep_mix.R')
source('main_estep_newkernel.R')
source('main_mamyula.R')
source('main_estep_mala.R')
source("mixtureFunctions.R")
library("mlxR")
library(sgd)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")


# Doc
# data(theo.saemix)
# theo.saemix_less <- theo.saemix[1:120,]
# # theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# saemix.data<-saemixData(name.data=theo.saemix_less,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")



library(saemix)
PD1.saemix<-read.table( "data/PD1.saemix.tab",header=T,na=".")
PD2.saemix<-read.table( "data/PD2.saemix.tab",header=T,na=".")
saemix.data1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))

PD2.saemix <- PD2.saemix[1:168,]
saemix.data2<-saemixData(name.data=PD2.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))



modelemax<-function(psi,id,xidep) {
# input:
# psi : matrix of parameters (3 columns, E0, Emax, EC50)
# id : vector of indices
# xidep : dependent variables (same nb of rows as length of id)
# returns:
# a vector of predictions of length equal to length of id
dose<-xidep[,1]
e0<-psi[id,1]
emax<-psi[id,2]
e50<-psi[id,3]
f<-e0+emax*dose/(e50+dose)
return(f)
}

saemix.model<-saemixModel(model=modelemax,description="Emax model",
psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,
dimnames=list(NULL,c("E0","Emax","EC50"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
byrow=TRUE),error.model="constant")


K1 = 100
K2 = 50
iterations = 1:(K1+K2+1)
gd_step = 0.01

end = K1+K2

seed0 = 39546
#RWM
options<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,0,0,0,0,0), nbiter.saemix = c(K1,K2))
theo_ref<-data.frame(saemix_mamyula(saemix.model,saemix.data2,options))
theo_ref <- cbind(iterations, theo_ref)

theo_ref[end,]

graphConvMC_twokernels(theo_ref,theo_ref, title="new kernel")
#saem with mala
options.mala<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5,0,0),nbiter.saemix = c(K1,K2),sigma.val = 0.01,gamma.val=0.01)
theo_mala<-data.frame(saemix_mamyula(saemix.model,saemix.data2,options.mala))
theo_mala <- cbind(iterations, theo_mala)


#saem with mamyula
options.mamyula<-list(seed=39546,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,0,6,0),nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01,lambda.val=0.2)
theo_mamyula<-data.frame(saemix_mamyula(saemix.model,saemix.data2,options.mamyula))
theo_mamyula <- cbind(iterations, theo_mamyula)

graphConvMC_twokernels(theo_ref,theo_mala, title="new kernel")
graphConvMC_twokernels(theo_ref,theo_mamyula, title="new kernel")
graphConvMC_threekernels(theo_ref,theo_mala,theo_mamyula, title="new kernel")


replicate = 15
seed0 = 395246

#RWM
final_rwm <- 0
for (j in 1:replicate){
  print(j)
  options<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,0,0,0,0,0), nbiter.saemix = c(K1,K2))
  theo_ref<-data.frame(saemix_mamyula(saemix.model,saemix.data2,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- j
  final_rwm <- rbind(final_rwm,theo_ref)
}


names(final_rwm)[1]<-paste("time")
names(final_rwm)[9]<-paste("id")
final_rwm1 <- final_rwm[c(9,1,2)]
final_rwm2 <- final_rwm[c(9,1,3)]
final_rwm3 <- final_rwm[c(9,1,4)]
final_rwm4 <- final_rwm[c(9,1,5)]
final_rwm5 <- final_rwm[c(9,1,6)]
final_rwm6 <- final_rwm[c(9,1,7)]
final_rwm7 <- final_rwm[c(9,1,8)]
# prctilemlx(final_rwm1[-1,],band = list(number = 8, level = 80)) + ggtitle("RWM")
final_mix <- 0
for (j in 1:replicate){
  print(j)
  options.mala<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5,0,0),nbiter.saemix = c(K1,K2),sigma.val = 0.01,gamma.val=0.01)
  theo_mix<-data.frame(saemix_mamyula(saemix.model,saemix.data2,options.mala))
  theo_mix <- cbind(iterations, theo_mix)
  theo_mix['individual'] <- j
  final_mix <- rbind(final_mix,theo_mix)
}

#mix (RWM and MAP new kernel for liste of saem iterations)
final_mix <- 0
for (j in 1:replicate){
  print(j)
  options.mamyula<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(0,0,0,0,6,0),nbiter.saemix = c(K1,K2),sigma.val = 0.1,gamma.val=0.01,lambda.val=0.2)
  theo_mix<-data.frame(saemix_mamyula(saemix.model,saemix.data2,options.mamyula))
  theo_mix <- cbind(iterations, theo_mix)
  theo_mix['individual'] <- j
  final_mix <- rbind(final_mix,theo_mix)
}


names(final_mix)[1]<-paste("time")
names(final_mix)[9]<-paste("id")
final_mix1 <- final_mix[c(9,1,2)]
final_mix2 <- final_mix[c(9,1,3)]
final_mix3 <- final_mix[c(9,1,4)]
final_mix4 <- final_mix[c(9,1,5)]
final_mix5 <- final_mix[c(9,1,6)]
final_mix6 <- final_mix[c(9,1,7)]
final_mix7 <- final_mix[c(9,1,8)]

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
                                    colors       = c('red', 'blue'))
plot.S <- plot.S1  + ylab("ka")+ theme(legend.position=c(0.9,0.8))+ theme_bw()
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
                                    colors       = c('red', 'blue'))
plot.S2 <- plot.S2  + ylab("V")+ theme(legend.position=c(0.9,0.8))+ theme_bw()


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
                                    colors       = c('red', 'blue'))
plot.S3 <- plot.S3  + ylab("k")+ theme(legend.position=c(0.9,0.8))+ theme_bw()




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
                                    colors       = c('red', 'blue'))
plot.S4 <- plot.S4  + ylab("w2ka")+ theme(legend.position=c(0.9,0.8))+ theme_bw()


final_rwm5['group'] <- 1
final_mix5['group'] <- 2
final_mix5$id <- final_mix5$id +1


final5 <- rbind(final_rwm5[-1,],final_mix5[-1,])
labels <- c("ref","new")
# prctilemlx(final5[c(1,4,2,3)], band = list(number = 4, level = 80),group='group', label = labels) 
# plt1 <- prctilemlx(final1, band = list(number = 4, level = 80),group='group', label = labels) 

# rownames(final1) <- 1:nrow(final1)

plot.S5 <- plot.prediction.intervals(final5[c(1,4,2,3)], 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('red', 'blue'))
plot.S5 <- plot.S5  + ylab("w2V")+ theme(legend.position=c(0.9,0.8))+ theme_bw()



final_rwm6['group'] <- 1
final_mix6['group'] <- 2
final_mix6$id <- final_mix6$id +1


final6 <- rbind(final_rwm6[-1,],final_mix6[-1,])
labels <- c("ref","new")
# prctilemlx(final6[c(1,4,2,3)], band = list(number = 4, level = 80),group='group', label = labels) 
# plt1 <- prctilemlx(final1, band = list(number = 4, level = 80),group='group', label = labels) 

# rownames(final1) <- 1:nrow(final1)

plot.S6 <- plot.prediction.intervals(final6[c(1,4,2,3)], 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('red', 'blue'))
plot.S6 <- plot.S6  + ylab("w2k")+ theme(legend.position=c(0.9,0.8))+ theme_bw()



final_rwm7['group'] <- 1
final_mix7['group'] <- 2
final_mix7$id <- final_mix7$id +1


final7 <- rbind(final_rwm7[-1,],final_mix7[-1,])
labels <- c("ref","new")
# prctilemlx(final7[c(1,4,2,3)], band = list(number = 4, level = 80),group='group', label = labels) 
# plt1 <- prctilemlx(final1, band = list(number = 4, level = 80),group='group', label = labels) 

# rownames(final1) <- 1:nrow(final1)

plot.S7 <- plot.prediction.intervals(final7[c(1,4,2,3)], 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('red', 'blue'))
plot.S7 <- plot.S7  + ylab("a")+ theme(legend.position=c(0.9,0.8))+ theme_bw()







grid.arrange(plot.S, plot.S2,plot.S3,plot.S4, plot.S5,plot.S6,plot.S7,ncol=3)


#values table

#values table
sample_mean_rwm <- 0
var_rwm <- 0
error_rwm <- 0
true_param <- c(1.5,32,0.1,0.4,0.01,0.8)
for (j in 1:replicate){
  sample_mean_rwm <- sample_mean_rwm + colMeans(final_rwm[(j*K1):(j*(K1+K2)),c(2,3,4,5,6,8)])
}
sample_mean_rwm = 1/replicate*sample_mean_rwm

for (j in 1:replicate){
  var_rwm <- var_rwm + (final_rwm[(j*(K1+K2)),c(2,3,4,5,6,8)]-sample_mean_rwm)^2
  error_rwm <- error_rwm + (final_rwm[(j*(K1+K2)),c(2,3,4,5,6,8)]-true_param)^2
}

error_rwm = 1/replicate*error_rwm
var_rwm = 1/replicate*var_rwm




sample_mean_mix <- 0
var_mix <- 0
error_mix <- 0
true_param <- c(1.5,32,0.1,0.4,0.01,0.8)
for (j in 1:replicate){
  sample_mean_mix <- sample_mean_mix + colMeans(final_mix[(j*K1):(j*(K1+K2)),c(2,3,4,5,6,8)])
}
sample_mean_mix = 1/replicate*sample_mean_mix

for (j in 1:replicate){
  var_mix <- var_mix + (final_mix[(j*(K1+K2)),c(2,3,4,5,6,8)]-sample_mean_mix)^2
  error_mix <- error_mix + (final_mix[(j*(K1+K2)),c(2,3,4,5,6,8)]-true_param)^2
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




