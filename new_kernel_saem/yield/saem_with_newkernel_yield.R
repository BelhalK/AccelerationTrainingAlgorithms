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
# source('initalgo.R')
source('main_new.R')
source('main_estep_new.R')
source('main_gd.R')
source('main_estep_gd.R')
source('main_gd_mix.R')
source('main_estep_gd_mix.R')
source('main_estep_mix.R')
source('main_estep_newkernel.R')
source("mixtureFunctions.R")
library("mlxR")
#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")


# Doc
data(yield.saemix)
saemix.data<-saemixData(name.data=yield.saemix,header=TRUE,name.group=c("site"),
  name.predictors=c("dose"),name.response=c("yield"),
  name.covariates=c("soil.nitrogen"),units=list(x="kg/ha",y="t/ha",
  covariates=c("kg/ha")))

yield.LP<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, ymax, xmax, slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  ymax<-psi[id,1]
  xmax<-psi[id,2]
  slope<-psi[id,3]
  f<-ymax+slope*(x-xmax)
#  cat(length(f),"  ",length(ymax),"\n")
  f[x>xmax]<-ymax[x>xmax]
  return(f)
}

saemix.model<-saemixModel(model=yield.LP,description="Linear plus plateau model",   
  psi0=matrix(c(8,100,0.2,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("Ymax","Xmax","slope"))),covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")



K1 = 70
K2 = 50
iterations = 1:(K1+K2+1)
gd_step = 0.00001
seed0 = 39546

# #RWM
options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0)
theo_ref<-data.frame(saemix_new(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)

# graphConvMC_twokernels(theo_ref,theo_ref, title="RWM vs Laplace SAEM")


#ref (map always)
options.new<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5),nbiter.saemix = c(K1,K2),nbiter.sa=0)
theo_new_ref<-data.frame(saemix_new(saemix.model,saemix.data,options.new))
theo_new_ref <- cbind(iterations, theo_new_ref)


graphConvMC_twokernels(theo_ref,theo_new_ref, title="RWM vs Laplace SAEM")


K1 = 100
K2 = 50
iterations = 1:(K1+K2+1)
gd_step = 0.01


# #mix (RWM and MAP new kernel for liste of saem iterations)
# options.mix<-list(seed=395246,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,4,0),nbiter.saemix = c(K1,K2),step.gd=gd_step,map.range=3)
# theo_mix<-data.frame(saemix_gd_mix(saemix.model,saemix.data,options.mix))
# theo_mix <- cbind(iterations, theo_mix)



replicate = 2
seed0 = 39546

#RWM
final_rwm <- 0
for (j in 1:replicate){
  print(j)
  options<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0)
  theo_ref<-data.frame(saemix_new(saemix.model,saemix.data,options))
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
final_rwm6 <- final_rwm[c(9,1,8)]
# prctilemlx(final_rwm1[-1,],band = list(number = 8, level = 80)) + ggtitle("RWM")

#mix (RWM and MAP new kernel for liste of saem iterations)
final_mix <- 0
for (j in 1:replicate){
  print(j)
  options.mix<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(1,0,0,5),nbiter.saemix = c(K1,K2),step.gd=gd_step,map.range=c(1:3))
  theo_mix<-data.frame(saemix_new(saemix.model,saemix.data,options.mix))
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
final_mix6 <- final_mix[c(9,1,8)]

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
plot.S <- plot.S1  + ylab("Ymax")+ theme(legend.position=c(0.9,0.8))+ theme_bw()
print(plot.S1)



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
plot.S2 <- plot.S2  + ylab("Xmax")+ theme(legend.position=c(0.9,0.8))+ theme_bw()


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
plot.S3 <- plot.S3  + ylab("slope")+ theme(legend.position=c(0.9,0.8))+ theme_bw()




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
plot.S4 <- plot.S4  + ylab("w1")+ theme(legend.position=c(0.9,0.8))+ theme_bw()


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
                                    colors       = c('#01b7a5', '#c17b01'))
plot.S5 <- plot.S5  + ylab("w2")+ theme(legend.position=c(0.9,0.8))+ theme_bw()



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
                                    colors       = c('#01b7a5', '#c17b01'))
plot.S6 <- plot.S6  + ylab("a")+ theme(legend.position=c(0.9,0.8))+ theme_bw()






grid.arrange(plot.S, plot.S2,plot.S3,plot.S4, plot.S5,plot.S6,ncol=3)

