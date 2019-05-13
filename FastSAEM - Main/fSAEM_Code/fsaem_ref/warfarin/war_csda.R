setwd("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/Dir")
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
  source('main_new.R')
  source('main_estep_new2.R')
  source('main_new_mix.R')
  source('main_estep_mix.R')
  
setwd("/Users/karimimohammedbelhal/Desktop/CSDA_code/")
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



setwd("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin")
# source('dataproc.R')
# model 
model<-"warfarin_project_model.txt"
# treatment
trt <- read.table("treatment.txt", header = TRUE) 

# parameters 
originalId<- read.table('originalId.txt', header=TRUE) 
populationParameter<- read.vector('populationParameter.txt') 
individualCovariate<- read.table('individualCovariate.txt', header = TRUE) 
list.param <- list(populationParameter,individualCovariate)
# output 
name<-"y1"
time<-read.table("output1.txt",header=TRUE)
out1<-list(name=name,time=time) 
name<-"y2"
time<-read.table("output2.txt",header=TRUE)
out2<-list(name=name,time=time) 
out<-list(out1,out2)

# call the simulator 
res <- simulx(model=model,treatment=trt,parameter=list.param,output=out)
# warfarin.saemix <- data(res)

# writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/warfarin/war_synth.csv")
# table <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/warfarin/war_synth.csv", header=T, sep=",")
# table <- table[table$ytype==1,]

# table[,5] <- 0

warfarin.saemix <- res$y1
warfarin.saemix["amount"] <- 0
treat <- res$treatment
treat["y1"] <- 0
treat <- treat[c(1,2,4,3)]

j <- 1
l<-c()


for (i in 1:241) {
    
    if(t(warfarin.saemix["id"])[i]==t(treat["id"])[j]){
        print(rownames(warfarin.saemix[i,]))
        l <- rbind(l,rownames(warfarin.saemix[i,]))
        j<-j+1
      }
}

warfarin.saemix <- rbind(treat[1,], warfarin.saemix)
warfarin.saemix[1:7,4] <- treat[1,4]
j <- 2
for (i in l[-1]){
  warfarin.saemix[(as.numeric(i)+1):(as.numeric(i)+length(which(t(warfarin.saemix["id"]==j)))),4] <- treat[j,4]
  warfarin.saemix <- rbind(warfarin.saemix[1:(as.numeric(i)-1),], treat[j,], warfarin.saemix[(as.numeric(i)+1):nrow(warfarin.saemix),])
  j <- j +1
}

rownames(warfarin.saemix) <- 1:nrow(warfarin.saemix)
# warfarin.saemix <- table[c(1,2,3,5)]
warfarin.saemix_less <- warfarin.saemix[,]

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  k<-psi[id,3]
  CL<-k*V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}

# warfa<-function(psi,id,xidep) { 
#   dose<-xidep[,1]
#   tim<-xidep[,2]  
#   ka<-psi[id,1]
#   V<-psi[id,2]
#   k<-psi[id,3]
#   CL<-k*V
#   Tlag<-psi[id,4]
#   Imax<-psi[id,5]
#   IC50<-psi[id,6]
#   kin<-psi[id,7]
#   kout<-psi[id,8]
#   ypred1<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#   E_0 <- kin/kout
#   dypred2<-kin*(1-Imax*ypred1/(ypred1+IC50))-kout*ypred2
#   return(ypred1,ypred2)
# }



# # Default model, no covariate
# saemix.model<-saemixModel(model=model1cpt,description="warfarin"
#   ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
#   transform.par=c(1,1,1))



# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="warfarin"
  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE))


saemix.data<-saemixData(name.data=warfarin.saemix_less,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")

K1 = 400
K2 = 200
iterations = 1:(K1+K2+1)
gd_step = 0.01
end = K1+K2
seed0 = 395246

#RWM
options<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,0,0),nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0)
theo_ref<-data.frame(saemix_new_mix(saemix.model,saemix.data,options))
theo_ref <- cbind(iterations, theo_ref)
theo_ref[end,]
graphConvMC_twokernels(theo_ref,theo_ref, title="new kernel")

#ref (map always)
# options.new<-list(seed=395246,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(0,0,0,6),nbiter.saemix = c(K1,K2))
# theo_new_ref<-data.frame(saemix_new(saemix.model,saemix.data,options.new))
# theo_new_ref <- cbind(iterations, theo_new_ref)



# theo_new_ref[end,]


options.new<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),map.range=c(1:10),nbiter.burn =0)
theo_new_ref<-data.frame(saemix_new_mix(saemix.model,saemix.data,options.new))
theo_new_ref <- cbind(iterations, theo_new_ref)

graphConvMC_twokernels(theo_ref,theo_new_ref, title="new kernel")
#mix (RWM and MAP new kernel for liste of saem iterations)
# options.mix<-list(seed=395246,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,4,0),nbiter.saemix = c(K1,K2),step.gd=gd_step,map.range=3)
# theo_mix<-data.frame(saemix_gd_mix(saemix.model,saemix.data,options.mix))
# theo_mix <- cbind(iterations, theo_mix)



#First run on the same dataset

replicate = 3

final_rwm <- 0
final_mix <- 0
for (m in 1:replicate){
  print(m)
  print(m)
  l = list(c(1,5,1,0,0,0),c(0.8,4,0.8,0,0,0),c(1.2,3,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  saemix.model<-saemixModel(model=model1cpt,description="warfarin"
  ,psi0=matrix(l[[m]],ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1/m,0,0,0,1/m,0,0,0,1/m),ncol=3,byrow=TRUE))

  options<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0)
  theo_ref<-data.frame(saemix_new_mix(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  theo_ref['individual'] <- m
  final_rwm <- rbind(final_rwm,theo_ref[-1,])

  options.new<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),map.range=c(1:10),nbiter.burn =0)
  theo_new_ref<-data.frame(saemix_new_mix(saemix.model,saemix.data,options.new))
  theo_mix <- cbind(iterations, theo_new_ref)
  theo_mix['individual'] <- m
  final_mix <- rbind(final_mix,theo_mix[-1,])
}

graphConvMC_diff2(final_rwm,final_mix, title="Diff intial param Warfa")
graphConvMC_diff2(final_rwm[,c(1,3,6,9)],final_mix[,c(1,3,6,9)], title="Diff intial param Warfa")

graphConvMC_diff(final_rwm[,c(1,2,9)],final_mix[,c(1,2,9)], title="Diff intial param Warfa")

graphConvMC_diff2(final_rwm[,c(1,3,6,9)],final_mix[,c(1,3,6,9)])



graphConvMC_diff3 <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=2) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=2)+
      xlab("") +scale_x_log10()+ ylab("")  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=14, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=14, angle=0))+ggtitle(names(df[j]))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


graphConvMC_diff4 <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=2) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=2)+
      xlab("") +scale_x_log10()+scale_y_continuous( limits=c(2, 8))+ ylab("")  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=14, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=14, angle=0))+ggtitle(names(df[j]))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


a <- graphConvMC_diff4(final_rwm[,c(1,3,9)],final_mix[,c(1,3,9)])
b <- graphConvMC_diff3(final_rwm[,c(1,6,9)],final_mix[,c(1,6,9)])

grid.arrange(a,b, ncol=2)

#run on diff datasets to compute residual errors

var_rwm <- 0
error_rwm <- 0


var_mix <- 0
error_mix <- 0

replicate = 4

final_rwm <- 0
final_mix <- 0

replicate <- 40
for (m in 1:replicate){
  
    model<-"warfarin_project_model.txt"
# treatment
trt <- read.table("treatment.txt", header = TRUE) 

# parameters 
originalId<- read.table('originalId.txt', header=TRUE) 
populationParameter<- read.vector('populationParameter.txt') 
individualCovariate<- read.table('individualCovariate.txt', header = TRUE) 
list.param <- list(populationParameter,individualCovariate)
# output 
name<-"y1"
time<-read.table("output1.txt",header=TRUE)
out1<-list(name=name,time=time) 
name<-"y2"
time<-read.table("output2.txt",header=TRUE)
out2<-list(name=name,time=time) 
out<-list(out1,out2)

# call the simulator 
res <- simulx(model=model,treatment=trt,parameter=list.param,output=out)
# warfarin.saemix <- data(res)

# writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/warfarin/war_synth.csv")
# table <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/warfarin/war_synth.csv", header=T, sep=",")
# table <- table[table$ytype==1,]

# table[,5] <- 0

warfarin.saemix <- res$y1
warfarin.saemix["amount"] <- 0
treat <- res$treatment
treat["y1"] <- 0
treat <- treat[c(1,2,4,3)]

j <- 1
l<-c()


for (i in 1:241) {
    
    if(t(warfarin.saemix["id"])[i]==t(treat["id"])[j]){
        print(rownames(warfarin.saemix[i,]))
        l <- rbind(l,rownames(warfarin.saemix[i,]))
        j<-j+1
      }
}

warfarin.saemix <- rbind(treat[1,], warfarin.saemix)
warfarin.saemix[1:7,4] <- treat[1,4]
j <- 2
for (i in l[-1]){
  warfarin.saemix[(as.numeric(i)+1):(as.numeric(i)+length(which(t(warfarin.saemix["id"]==j)))),4] <- treat[j,4]
  warfarin.saemix <- rbind(warfarin.saemix[1:(as.numeric(i)-1),], treat[j,], warfarin.saemix[(as.numeric(i)+1):nrow(warfarin.saemix),])
  j <- j +1
}

rownames(warfarin.saemix) <- 1:nrow(warfarin.saemix)
# warfarin.saemix <- table[c(1,2,3,5)]
warfarin.saemix_less <- warfarin.saemix[,]
  saemix.model<-saemixModel(model=model1cpt,description="warfarin"
  ,psi0=matrix(c(1,4,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
  transform.par=c(1,1,1),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE))


saemix.data<-saemixData(name.data=warfarin.saemix_less,header=TRUE,sep=" ",na=NA, name.group=c("id"),
  name.predictors=c("amount","time"),name.response=c("y1"), name.X="time")




  options<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),map.range=c(0),nbiter.burn =0)
  theo_ref<-data.frame(saemix_new_mix(saemix.model,saemix.data,options))
  theo_ref <- cbind(iterations, theo_ref)
  # var_rwm <- var_rwm + (theo_ref[,2:8]-true_param)^2
  ML <- theo_ref[,2:8]
  ML[1:(end+1),]<- theo_ref[end+1,2:8]
  error_rwm <- error_rwm + (theo_ref[,2:8]-ML)^2
  theo_ref['individual'] <- j
  final_rwm <- rbind(final_rwm,theo_ref)
  

  options.new<-list(seed=seed0,map=F,fim=F,ll.is=T,nb.chains = 1, nbiter.mcmc = c(2,2,2,6,0),nbiter.saemix = c(K1,K2),map.range=c(1:3),nbiter.burn =0)
  print(m)
  theo_mix<-data.frame(saemix_new_mix(saemix.model,saemix.data,options.new))
  theo_mix <- cbind(iterations, theo_mix)
  # var_mix <- var_mix + (theo_mix[,2:8]-true_param)^2
  ML <- theo_mix[,2:8]
  ML[1:(end+1),]<- theo_mix[end+1,2:8]
  error_mix <- error_mix + (theo_mix[,2:8]-ML)^2
  theo_mix['individual'] <- m
  final_mix <- rbind(final_mix,theo_mix)
}

error_rwm <- 1/replicate*error_rwm
error_mix <- 1/replicate*error_mix

error_rwm <- cbind(iterations, error_rwm)
error_mix <- cbind(iterations, error_mix)

err_mix<- theo_ref
err_rwm<- theo_ref
err_rwm[,2:8] <- error_rwm[,2:8]
err_mix[,2:8] <- error_mix[,2:8]

err_mix[2,] = err_rwm[2,]


# graphConvMC_diff(err_rwm[-1,],err_mix[-1,], title="Quadratic errors Warfa")
# graphConvMC_diff2(err_rwm[-1,],err_mix[-1,], title="Quadratic errors Warfa")

# graphConvMC_diff2(err_rwm[-1,c(1,3,6,9)],err_mix[-1,c(1,3,6,9)])

a <- err_rwm[-1,]
b <- err_mix[-1,]
# graphConvMC_diff2(err_rwm[-1,c(1,3,6,9)],err_mix[-1,c(1,3,6,9)], title="Quadratic errors Warfa")

# err_rwm[,2:8] <- sqrt(err_rwm[,2:8])
# err_mix[,2:8] <- sqrt(err_mix[,2:8])

# graphConvMC_diff2(a[,c(1,3,6,9)],b[,c(1,3,6,9)], title="Quadratic errors Warfa")
# graphConvMC_diff2(err_rwm[-1,c(1,3,6,9)],err_mix[-1,c(1,3,6,9)], title="Quadratic errors Warfa")

# graphConvMC_diff2(a[,c(1,3,6,9)],b[,c(1,3,6,9)])

c <- graphConvMC_diff3(a[,c(1,3,9)],b[,c(1,3,9)])
d <- graphConvMC_diff3(a[,c(1,6,9)],b[,c(1,6,9)])

grid.arrange(c,d, ncol=2)

library(FactoMineR)
data(decathlon)
X <- decathlon[, 1:10]
res_pca <- PCA(X, graph = F)
ind_summary <- res_pca$ind$coord[,1:2]
head(ind_summary)

