source("algos.R")
source("func.R")
theme_set(theme_bw())
# save.image("lin_gauss.RData")
load("lin_gauss.RData")
# n <- 100
# mu<-c(1,0)
# mu0<-c(0.1,0)
# sigma<-c(0.5,0.1)*1

n <- 300
mu<-c(10,0)
mu0<-c(5,0)
# sigma<-c(10,5)*1
sigma<-c(1,1)*1

alph <- sigma[2]/(sigma[1]+sigma[2])
gamm <- 1/(1/sigma[1]+1/sigma[2])

K <- 3000
# Several Chains for the same iteration
M <- 1

# alpha1 <- 0.7
# alpha2 <- 0.4
seed0=44444
ylim <- c(0.3)

nsim <- 2
G<-1
col.names <- c("iteration", paste0("mu",1:G))
theta<-list(mu=mu[1])
theta0<-list(mu=mu0[1])


set.seed(seed0)
xj<-mixt.simulate(n,mu,sigma)
df <- mixt.em(xj, theta0, K, alph)

ML <- df
for (i in (1:(K+1))){
  ML[i,2]<- c(theta$mu)
}

## EM
print('EM')
dem <- NULL
nbrem<-n
df.em <- vector("list", length=nsim)

nbr<-1
diem <- NULL
df.iem <- vector("list", length=nsim)

doem <- NULL
df.oem <- vector("list", length=nsim)

doemvr <- NULL
df.oemvr <- vector("list", length=nsim)

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  x <- NULL
  for (i in 1:n){
   xj<-mixt.simulate(n,mu,sigma)
   x <- rbind(x,xj)
  }

  df <- mixt.em(x[,j], theta0, K, alph)
  # ML <- df
  # ML[1:(K+1),2]<- df[(K+1),2]
  df[,2] <- (df[,2] - ML[,2])^2
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df

  df <- mixt.iem(x[,j], theta0, K, alph,nbr)
  df[,2] <- (df[,2] - ML[,2])^2
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df

  df <- mixt.oem(x[,j], theta0, K, alph,nbr)
  df[,2] <- (df[,2] - ML[,2])^2
  df$rep <- j
  doem <- rbind(doem,df)
  df$rep <- NULL
  df.oem[[j]] <- df

  df <- mixt.oemvr(x[,j], theta0, K, alph,nbr)
  df[,2] <- (df[,2] - ML[,2])^2
  df$rep <- j
  doemvr <- rbind(doemvr,df)
  df$rep <- NULL
  df.oemvr[[j]] <- df
}


# dem[,2] <- dem[,2]^2
em <- NULL
em <- dem[dem$rep==1,]

if (nsim>2) {
   for (j in (2:nsim))
	{
	  em[,2] <- em[,2]+dem[dem$rep==j,2]
	}
}
em[,2] <- 1/nsim*em[,2]
em[,4]<-NULL



iem <- NULL
iem <- diem[diem$rep==1,]

if (nsim>2) {
		for (j in (2:nsim))
	{
	  iem[,2] <- iem[,2]+diem[diem$rep==j,2]
	}
}

iem[,2] <- 1/nsim*iem[,2]
iem[,4]<-NULL



oem <- NULL
oem <- doem[doem$rep==1,]

if (nsim>2) {
    for (j in (2:nsim))
  {
    oem[,2] <- oem[,2]+doem[doem$rep==j,2]
  }
}

oem[,2] <- 1/nsim*oem[,2]
oem[,4]<-NULL



oemvr <- NULL
oemvr <- doemvr[doemvr$rep==1,]

if (nsim>2) {
    for (j in (2:nsim))
  {
    oemvr[,2] <- oemvr[,2]+doemvr[doemvr$rep==j,2]
  }
}

oemvr[,2] <- 1/nsim*oemvr[,2]
oemvr[,4]<-NULL



em_scaled <- em
em_scaled$iteration = seq(0, n*K, by=n)
em_scaled <- em_scaled[rep(seq_len(nrow(em_scaled)), each=n),]


iem$algo <- 'IEM'
oem$algo <- 'OEM'
oemvr$algo <- 'OEMvr'
em_scaled$algo <- 'EM'
# variance <- NULL
# variance <- rbind(em_scaled[0:K,],iem[0:K,],oem[0:K,],oemvr[0:K,])
# colnames(variance) <- c("iteration","mu1","algo")

variance <- NULL
variance <- rbind(em_scaled[0:(K+1),c(1,2,4)],iem[0:(K+1),c(1,2,4)],oem[0:(K+1),c(1,2,4)],oemvr[0:(K+1),c(1,2,4)])
graphConvMC2_new(variance, title="IEMs",legend=TRUE)



# graphConvMC <- function(df,df2,df3,df4, title=NULL, ylim=NULL)
# {
#   G <- (ncol(df)-2)/3
#   df$algo <- as.factor(df$algo)
#   df2$algo <- as.factor(df2$algo)
#   df3$algo <- as.factor(df3$algo)
#   df4$algo <- as.factor(df4$algo)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",linetype= "solid",size=2)+
#     geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="black",linetype="longdash",size=2)+
#     geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red",linetype="dotted",size=2)+
#     geom_line(aes_string(df4[,1],df4[,j],by=df4[,ncol(df4)]),colour="blue",linetype="dotted",size=2)+
#       xlab("") +ylab(expression(paste(beta,"1")))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
#                            size=30, angle=0),
#           axis.text.y = element_text(face="bold", color="black", 
#                            size=30, angle=0))+theme(axis.title = element_text(color="black", face="bold", size=30)) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj

#   }
#   do.call("grid.arrange", c(graf, ncol=1, top=title))
# }

# m <- graphConvMC(em_scaled[0:K,c(1,2,4)],iem[0:K,c(1,2,4)],oem[0:K,c(1,2,4)],oemvr[0:K,c(1,2,4)])


