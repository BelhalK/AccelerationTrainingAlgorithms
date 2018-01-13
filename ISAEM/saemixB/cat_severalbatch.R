setwd("/Users/karimimohammedbelhal/Desktop/CSDA_code/Dir")
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
  source('main_time.R')
  source('main_estep_time.R')
  source('main_mstep_time.R') 
  source('func_aux_time.R') 
  source('SaemixObject_time.R') 
  source('main_initialiseMainAlgo_time.R') 
setwd("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat_incremental")
source('main_cat2.R')
source('main_estep_cat2.R')
source('main_cat_incremental.R')
source('estep_cat_incremental.R')
source('main_mstep_cat.R') 
source('func_aux_cat.R') 
source('SaemixObject_cat.R') 
source('main_initialiseMainAlgo_cat.R') 
source("mixtureFunctions.R")


library("rJava")
library("rCMA")
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
iter_mcmc = 200


# cat_data.saemix<-read.table("data/categorical1_data.txt",header=T,na=".")
# cat_data.saemix <- cat_data.saemix[1:692,]
# saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("ID"),name.response=c("Y"),name.predictors=c("Y"), name.X=c("TIME"))



# cat_data.saemix<-read.table("data/categorical1_data.txt",header=T,na=".")
# cat_data.saemix<-read.table("data/categorical1_data_less2.txt",header=T,na=".")
# cat_data.saemix<-read.table("data/categorical1_data_less2.txt",header=T,na=".")
# cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat1.csv", header=T, sep=",")
cat_data.saemix<-read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat2.csv", header=T, sep=",")
saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("id"),name.response=c("y"),name.predictors=c("y"), name.X=c("time"))
# saemix.data<-saemixData(name.data=cat_data.saemix,header=TRUE,sep=" ",na=NA, name.group=c("ID"),name.response=c("Y"),name.predictors=c("Y"), name.X=c("TIME"))



cat_data.model<-function(psi,id,xidep) {
level<-xidep[,1]

th1 <- psi[id,1]
th2 <- psi[id,2]
th3 <- psi[id,3]

P0 <- 1/(1+exp(-th1))
Pcum1 <- 1/(1+exp(-th1-th2))
Pcum2 <- 1/(1+exp(-th1-th2-th3))

P1 <- Pcum1 - P0
P2 <- Pcum2 - Pcum1
P3 <- 1 - Pcum2

P.obs = (level==0)*P0+(level==1)*P1+(level==2)*P2+(level==3)*P3

return(P.obs)
}


saemix.model<-saemixModel(model=cat_data.model,description="cat model",   
  psi0=matrix(c(2,1,2),ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))), 
  transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),omega.init=matrix(c(2,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")


saemix.options_rwm<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(iter_mcmc,0,0,0))
saemix.foce<-list(seed=39546,map=F,fim=F,ll.is=F, nb.chains = 1, nbiter.mcmc = c(1,0,0,iter_mcmc))


# post_rwm<-saemix_post_cat(saemix.model,saemix.data,saemix.options_rwm)$post_rwm
# post_foce<-saemix_post_cat(saemix.model,saemix.data,saemix.foce)$post_newkernel


K1 = 500
K2 = 100

iteration = 1:(K1+K2+1)
gd_step = 0.01
end = K1+K2
seed0 = 444

#RWM
theo_ref <- NULL
options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=100)
theo_ref<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options))
theo_ref <- cbind(iteration, theo_ref)


# graphConvMC_saem(theo_ref, title="new kernel")


#RWM
cat_incremental <- NULL
options.incremental<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=50)
cat_incremental<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental))
cat_incremental <- cbind(iteration, cat_incremental)

cat_incremental25 <- NULL
options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=25)
cat_incremental25<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental25))
cat_incremental25 <- cbind(iteration, cat_incremental25)

graphConvMC2_saem(theo_ref,cat_incremental25, title="new kernel")

cat_incremental[end,]


theo_ref$algo <- 'MCEM'
cat_incremental$algo <- 'IMCEM50'
cat_incremental25$algo <- 'IMCEM25'


theo_ref_scaled <- theo_ref[rep(seq_len(nrow(theo_ref)), each=4),]
cat_incremental_scaled <- cat_incremental[rep(seq_len(nrow(cat_incremental)), each=2),]
# theo_ref_scaled <- theo_ref[rep(seq_len(nrow(theo_ref)), each=2),]
theo_ref_scaled$iteration = 1:(4*(K1+K2+1))
cat_incremental_scaled$iteration = 1:(2*(K1+K2+1))


comparison <- 0
# comparison <- rbind(theo_ref,theo_incremental)
# comparison <- rbind(theo_ref_scaled[iteration,],cat_incremental)
comparison <- rbind(theo_ref_scaled[iteration,],cat_incremental_scaled[iteration,], cat_incremental25)

var <- melt(comparison, id.var = c('iteration','algo'), na.rm = TRUE)
graphConvMC3_new(var, title="ALGO - EM (same complexity)",legend=TRUE)



replicate = 3

final_rwm <- 0
final_incremental <- 0
final_incremental25 <- 0
for (m in 1:replicate){
  print(m)
  print(m)
  # l = list(c(1,5,1,0,0,0),c(0.8,12,0.8,0,0,0),c(1.2,3,1.2,0,0,0),c(1.4,6.6,1.4,0,0,0))
  l = list(c(2,1,2),c(2,1,5),c(2,1,3))
  saemix.model<-saemixModel(model=cat_data.model,description="cat model",   
  psi0=matrix(l[[m]],ncol=3,byrow=TRUE,dimnames=list(NULL,   
  c("th1","th2","th3"))), 
  transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),omega.init=matrix(c(3/m,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")


  # options<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=100)
  # theo_ref<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options))
  # theo_ref <- cbind(iteration, theo_ref)
  # theo_ref['individual'] <- m
  # theo_ref_scaled <- theo_ref[rep(seq_len(nrow(theo_ref)), each=4),]
  # theo_ref_scaled$iteration = 1:(4*(K1+K2+1))
  # final_rwm <- rbind(final_rwm,theo_ref_scaled[iteration,])

  options.incremental<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=50)
  theo_mix<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental))
  theo_mix <- cbind(iteration, theo_mix)
  theo_mix['individual'] <- m
  theo_mix_scaled <- theo_mix[rep(seq_len(nrow(theo_mix)), each=2),]
  theo_mix_scaled$iteration = 1:(2*(K1+K2+1))
  final_incremental <- rbind(final_incremental,theo_mix_scaled[iteration,])

  #  options.incremental25<-list(seed=seed0,map=F,fim=F,ll.is=F,nb.chains = 1, nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),displayProgress=TRUE, map.range=c(0),nbiter.sa=0,nbiter.burn =0, nb.replacement=25)
  # theo_mix25<-data.frame(saemix_cat_incremental(saemix.model,saemix.data,options.incremental25))
  # theo_mix25 <- cbind(iteration, theo_mix25)
  # theo_mix25['individual'] <- m
  # final_incremental25 <- rbind(final_incremental25,theo_mix25[,])
}

graphConvMC_diff3(final_incremental,final_incremental, final_incremental,title="Diff intial param Warfa")

graphConvMC_diff2(final_rwm[,c(1,3,9)],final_incremental[,c(1,3,9)], title="Diff intial param Warfa")

graphConvMC_diff2(final_rwm[,c(1,3,6,9)],final_incremental[,c(1,3,6,9)], title="Diff intial param Warfa")

graphConvMC_diff(final_rwm[,c(1,2,9)],final_incremental[,c(1,2,9)], title="Diff intial param Warfa")

graphConvMC_diff2(final_rwm[,c(1,3,6,9)],final_incremental[,c(1,3,6,9)])



graphConvMC_diff3 <- function(df,df2,df3, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  df3$individual <- as.factor(df3$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=2) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=2)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="green",linetype = 2,size=2)+
      xlab("") +scale_x_log10()+ ylab(expression(paste(omega,"3")))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=20, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=20, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=30)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}



graphConvMC_diff4 <- function(df,df2,df3, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$individual <- as.factor(df$individual)
  df2$individual <- as.factor(df2$individual)
  df3$individual <- as.factor(df3$individual)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=2) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=2)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="green",linetype = 2,size=2)+
      xlab("") +scale_x_log10()+ ylab(names(df[j]))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=20, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=20, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=30)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}


# a <- graphConvMC_diff4(final_rwm[,c(1,4,8)],final_incremental[,c(1,4,8)])
# b <- graphConvMC_diff3(final_rwm[,c(1,7,8)],final_incremental[,c(1,7,8)])

# grid.arrange(a,b, ncol=2)
colnames(final_rwm)[4] <- "z3"
colnames(final_incremental)[4] <- "z3"
colnames(final_incremental25)[4] <- "z3"
a <- graphConvMC_diff4(final_rwm[,c(1,4,8)],final_incremental[,c(1,4,8)],final_incremental25[,c(1,4,8)])
b <- graphConvMC_diff3(final_rwm[,c(1,5,8)],final_incremental[,c(1,5,8)],final_incremental25[,c(1,5,8)])

grid.arrange(a,b, ncol=2)

