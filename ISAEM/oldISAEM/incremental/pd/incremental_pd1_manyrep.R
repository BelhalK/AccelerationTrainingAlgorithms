
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
  
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/incremental")
source('main_incremental.R')
source('main_estep_incremental.R')
source("mixtureFunctions.R")


setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/pd")

library("mlxR")
library("psych")
library("coda")
library("Matrix")
#####################################################################################
# Theophylline

# Data - changing gender to M/F
# theo.saemix<-read.table("data/theo.saemix.tab",header=T,na=".")
# theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
# library(saemix)
PD1.saemix<-read.table( "PD1.saemix.tab",header=T,na=".")
PD2.saemix<-read.table( "PD2.saemix.tab",header=T,na=".")
saemix.data1<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"),
units=list(x="mg",y="-",covariates="-"))
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

K1 = 200
K2 = 30
iteration = 1:(K1+K2+1)
replicate = 2
seed0 = 39546

final_ref <- 0
for (j in 1:replicate){
  print("ref")
  print(j)
  options.ref<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2), nbiter.saemix = c(K1,K2),nbiter.sa=0,displayProgress=FALSE)
theo_ref<-data.frame(saemix(saemix.model,saemix.data1,options.ref))
theo_ref <- cbind(iteration, theo_ref)
theo_ref_scaled <- theo_ref[rep(seq_len(nrow(theo_ref)), each=2),]
theo_ref_scaled$iteration = 1:(2*(K1+K2+1))
theo_ref <- theo_ref_scaled[iteration,]

  theo_ref['individual'] <- j
  final_ref <- rbind(final_ref,theo_ref)
}



names(final_ref)[1]<-paste("time")
names(final_ref)[9]<-paste("id")
final_ref1 <- final_ref[c(9,1,4)]
# prctilemlx(final_mala1[-1,],band = list(number = 2, level = 80)) + ggtitle("mala")

#map always 
final_incremental <- 0
for (j in 1:replicate){
  print("incremental")
  print(j)
  options.incremental<-list(seed=j*seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2),nbiter.saemix = c(K1,K2),nb.replacement=50,displayProgress=FALSE)
  theo_incremental<-data.frame(incremental_saemix(saemix.model,saemix.data1,options.incremental))
  theo_incremental <- cbind(iteration, theo_incremental)
  theo_incremental['individual'] <- j
  final_incremental <- rbind(final_incremental,theo_incremental)
}


names(final_incremental)[1]<-paste("time")
names(final_incremental)[9]<-paste("id")
final_incremental1 <- final_incremental[c(9,1,4)]




final_ref1['group'] <- as.integer(1)
final_incremental1['group'] <- as.integer(2)
final_incremental1$id <- final_incremental1$id +1



final <- 0
final <- rbind(final_ref1[-1,],final_incremental1[-1,])



labels <- c("SAEM","ISAEM 50%")
final <- final[c(1,4,2,3)]
prctilemlx(final, band = list(number = 2, level = 80),group='group', label = labels, facet=FALSE) + theme(legend.position = "none")+ ggtitle(colnames(final)[4])


final$group <- as.integer(final$group)


plot.S1 <- plot.prediction.intervals(final, 
                                    labels       = labels, 
                                    legend.title = "Algorithm")
plot.S <- plot.S1  + ylab("EC50") + theme(legend.position=c(0.9,0.8)) + theme_bw()
print(plot.S1+ theme_bw())


plot.S1 <- plot.prediction.intervals(res1$S, 
                                    labels       = c("placebo","25 mg","50mg"), 
                                    legend.title = "arm")
plot.S <- plot.S1  + ylab("Survival prediction interval") + theme(legend.position=c(0.9,0.8))
print(plot.S1 )



library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)





plt <- prctilemlx(final, band = list(number = 4, level = 80),group='group', label = labels)
blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob
grid.arrange(plt,plot.S, ncol=2)


plot.prediction.intervals <- function(r, plot.median=TRUE, level=80, labels=NULL, 
                                      legend.title=NULL, colors=NULL) {
  P <- prctilemlx(r, number=1, level=level, plot=FALSE)
  if (is.null(labels))  labels <- a
  if (is.null(legend.title))  legend.title <- "group"
  names(P$y)[2:4] <- c("p.min","p50","p.max")
  pp <- ggplot(data=P$y)+ylab(NULL)+ 
    geom_ribbon(aes(x=time,ymin=p.min, ymax=p.max,fill=group),alpha=.5) 
  if (plot.median)
    pp <- pp + geom_line(aes(x=time,y=p50,colour=group))
  
  if (is.null(colors)) {
    pp <- pp + scale_fill_discrete(name=legend.title,
                                   breaks=a,
                                   labels=labels)
    pp <- pp + scale_colour_discrete(name=legend.title,
                                     breaks=a,
                                     labels=labels, 
                                     guide=FALSE)
  } else {
    pp <- pp + scale_fill_manual(name=legend.title,
                                 breaks=,
                                 labels=labels,
                                 values=colors)
    pp <- pp + scale_colour_manual(name=legend.title,
                                   breaks=a,
                                   labels=labels,
                                   guide=FALSE,values=colors)
  }  
  return(pp)
}

