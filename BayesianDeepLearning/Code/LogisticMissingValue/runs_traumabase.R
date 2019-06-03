library(MASS)
library(ggplot2)
library(reshape2)

source("miss.saem.R")
source("louis_lr_saem.R")
source("log_reg.R")
source("likelihood_saem.R")

# import
traumdata <- read.csv(file="trauma.csv",sep=';',dec=',',header = TRUE,na.strings = c("","NR","IMP","NA","NF"),encoding = "UTF-8")
SAMU <- traumdata[ ,c(10:14,225:228,234,33:35,15,18,229:233,244,19:30,41:42,49:54,48,58,61,44,46:47,55:57,60,62, 38)]

# preprocessing
SAMU$BMI[3302]=SAMU$Poids[3302]/(SAMU$Taille[3302]^2)

SAMU=SAMU[-which(SAMU$ACR.1==1),]
Choc.hemorragique= traumdata$Choc.hemorragique
Choc.hemorragique = Choc.hemorragique[-which(traumdata$ACR.1==1)]

Choc.hemorragique = Choc.hemorragique[-which(SAMU$Mecanisme=="Arme blanche" | SAMU$Mecanisme=="Arme à feu")]
SAMU=SAMU[-which(SAMU$Mecanisme=="Arme blanche" | SAMU$Mecanisme=="Arme à feu"),]

SAMU_CONT = SAMU[,c(1,3:5,34:41,43:44,48:50)]
indx <- sapply(SAMU_CONT, is.factor)
SAMU_CONT[indx] <- lapply(SAMU_CONT[indx], function(x) as.numeric(as.character(x)))
SAMU_CONT$SD.min=SAMU_CONT$PAS.min-SAMU_CONT$PAD.min
SAMU_CONT$SD.SMUR=SAMU_CONT$PAS.SMUR-SAMU_CONT$PAD.SMUR
SAMU_NEW = SAMU_CONT[,-c(7:8,10:11,15)]

nb.samples <- 600
nb.cov <- 10
# a subset of data to try
X.obs <- SAMU_NEW[1:nb.samples,1:nb.cov]
y <- Choc.hemorragique[1:nb.samples]


N = nrow(X.obs)
batchsize = 1
nb.epochs <-10
# nb.iter <- N/batchsize*nb.epochs
nb.iter <- 100

nb.epochs <-10
# SAEM
list.saem = miss.saem(X.obs,y,maxruns=10*nb.epochs,k1=0,ll_obs_cal=FALSE, algo = "saem")
#MCEM
list.mcem = miss.saem(X.obs,y,maxruns=10*nb.epochs,ll_obs_cal=FALSE, algo = "mcem")
#Incremental MCEM
list.imcem = miss.saem(X.obs,y,maxruns=N/batchsize*nb.epochs,ll_obs_cal=FALSE, algo = "imcem", batchsize= batchsize)
list.imcem10 = miss.saem(X.obs,y,maxruns=10*nb.epochs,ll_obs_cal=FALSE, algo = "imcem", batchsize= N/10)
list.imcem50 = miss.saem(X.obs,y,maxruns=2*nb.epochs,ll_obs_cal=FALSE, algo = "imcem", batchsize= N/2)

save(N, nb.iter, batchsize,list.saem,list.mcem,list.imcem,list.imcem10,list.imcem50, file = "trauma.RData")
