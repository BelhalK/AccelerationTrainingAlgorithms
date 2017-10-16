source("algos.R")
source("func.R")
theme_set(theme_bw())
library("mlxR")
# n <- 100
# mu<-c(1,0)
# mu0<-c(0.1,0)
# sigma<-c(0.5,0.1)*1

n <- 100
mu<-c(10,0)
mu0<-c(5,0)
# sigma<-c(10,5)*1
sigma<-c(3,10)*1

alph <- sigma[2]/(sigma[1]+sigma[2])
gamm <- 1/(1/sigma[1]+1/sigma[2])

KNR <- 300
# KR <- 400
K1<-30
# Several Chains for the same iteration
M <- 1

alpha1 <- 0.7
alpha2 <- 0.4
seed0=44444
ylim <- c(0.3)

nsim <- 5
G<-1
col.names <- c("iteration", paste0("mu",1:G))
theta<-list(mu=mu[1])
theta0<-list(mu=mu0[1])
##  Simulation
x <- matrix(0,nrow=n,ncol=nsim)
for (j in (1:nsim))
{
  seed <- j*seed0
  set.seed(seed)
  xj<-mixt.simulate(n,mu,sigma)
  x[,j] <- xj
}

## EM
print('EM')
dem <- NULL
nbr<-n
KR <- KNR*n/nbr
K<-KR
df.em <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.em(x[,j], theta0, KR, alph)
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df
}
# graphConvMC_new(dem, title="EM")


##########################################################################################


## MCEM
print('MCEM')
dmcem <- NULL
nbr<-n
KR <- KNR*n/nbr

df.mcem <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.mcem(x[,j], theta0, KR,K1,M=c(10,10),alph=alph,gamm=gamm)
  # for (k in (0:KR))
  # {
  #   df[k+1,2] <- (df[k+1,2] - df.em[[j]][1000,2])^2
  # }
  df$rep <- j

  df <- df[rep(seq_len(nrow(df)), each=10),]
  df$iteration = 1:(10*(KR+1))

  dmcem <- rbind(dmcem,df[1:KR,])
  df$rep <- NULL
  df.mcem[[j]] <- df
}




## IMCEM
print('IMCEM')
dimcem10 <- NULL
nbr <- n/10
KR <- KNR*n/nbr

df.imcem10 <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.imcem(x[,j], theta0, KR,K1, M=c(10,10),alph,gamm,nbr)
  # for (k in (0:KR))
  # {
  #   df[k+1,2] <- (df[k+1,2] - df.em[[j]][1000,2])^2
  # }
  df$rep <- j
  dimcem10 <- rbind(dimcem10,df)
  df$rep <- NULL
  df.imcem10[[j]] <- df
}



## IMCEM
print('IMCEM 50')
diffmcem50 <- NULL
nbr <- n/2
KR <- KNR*n/nbr

df.imcem50 <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  df <- mixt.imcem(x[,j], theta0, KR,K1, M=c(10,10),alph,gamm,nbr)
  # for (k in (0:KR))
  # {
  #   df[k+1,2] <- (df[k+1,2] - df.em[[j]][1000,2])^2
  # }
  df$rep <- j
  df <- df[rep(seq_len(nrow(df)), each=5),]
  df$iteration = 1:(5*(KR+1))

  diffmcem50 <- rbind(diffmcem50,df[1:KR,])
  df$rep <- NULL
  df.imcem50[[j]] <- df
}



dmcem['group'] <- 1
dmcem <- dmcem[c(3,4,1,2)]
colnames(dmcem)[3]<-"time"
colnames(dmcem)[1]<-"id"
diffmcem50['group'] <- 2
diffmcem50 <- diffmcem50[c(3,4,1,2)]
colnames(diffmcem50)[3]<-"time"
colnames(diffmcem50)[1]<-"id"
diffmcem50$id <- diffmcem50$id +1

dimcem10['group'] <- 3
dimcem10 <- dimcem10[c(3,4,1,2)]
colnames(dimcem10)[3]<-"time"
colnames(dimcem10)[1]<-"id"
dimcem10$id <- dimcem10$id +2


final <- rbind(dmcem,diffmcem50,dimcem10)

labels <- c("mcem","imcem 50%","imcem 10%")

prctilemlx(final, band = list(number = 4, level = 80),group='group', label = labels) 

plt <- prctilemlx(final, band = list(number = 4, level = 80),group='group', label = labels) 

rownames(final) <- 1:nrow(final)

plot.S1 <- plot.prediction.intervals(final, 
                                    labels       = labels, 
                                    legend.title = "algos",
                                    colors       = c('#01b7a5', '#c17b01', '#a00159'))
plot.S <- plot.S1  + ylab("mu1")+ theme(legend.position=c(0.9,0.8))+ theme_bw()
print(plot.S1)

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

grid.arrange(plt, plot.S, ncol=2)
# imcem$algo <- '10R'
imcem$algo <- '1R'
imcem50$algo <- '50R'
mcem$algo <- 'NR'
variance <- NULL
variance <- rbind(mcem[1000:30000,],imcem[1000:30000,], imcem50[1000:30000,])
variance <- rbind(mcem,imcem, imcem50)
graphConvMC2_new(variance, title="IMCEMs seq picking",legend=TRUE)


dmcem['group'] <- 1
final<-0
final <- dmcem[c(3,4,1,2)]
colnames(final)[3]<-"time"
colnames(final)[2]<-"id"
library("mlxR")
prctilemlx(final,band = list(number = 2, level = 80)) + ggtitle("mala")



dimcem['group'] <- 1
final<-0
final <- dimcem[c(3,4,1,2)]
colnames(final)[3]<-"time"
colnames(final)[2]<-"id"
library("mlxR")
prctilemlx(dimcem10,band = list(number = 2, level = 80)) + ggtitle("mala")



diffmcem50['group'] <- 1
final<-0
final <- diffmcem50[c(3,4,1,2)]
colnames(final)[3]<-"time"
colnames(final)[2]<-"id"
library("mlxR")
prctilemlx(final,band = list(number = 2, level = 80)) + ggtitle("mala")


# #####
# ## IMCEM
# print('IMCEM')
# dimcem <- NULL
# nbr <- n/100
# ind1 <- NULL
# KR <- KNR*n/nbr
# df.imcem <- vector("list", length=nsim)
# for (j in (1:nsim))
# {
#   print(j)
#   df <- mixt.imcem(x[,j], theta0, KR,K1, M=c(10,10),alph,gamm,nbr)
#   ind1 <- rbind(ind1,df$indices)
#   df<-df$param
#   # for (k in (0:KR))
#   # {
#   #   df[k+1,2] <- (df[k+1,2] - df.em[[j]][1000,2])^2
#   # }
#   df$rep <- j
#   dimcem <- rbind(dimcem,df)
#   df$rep <- NULL
#   df.imcem[[j]] <- df
# }

# ind1 <- as.data.frame(table(ind1))
# ind1.freq= as.vector(rep(ind1$ind1, ind1$Freq))
# indices3 <- as.numeric(ind1.freq)
# histo2 <- hist(indices3)


plot.prediction.intervals <- function(r, plot.median=TRUE, level=90, labels=NULL, 
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