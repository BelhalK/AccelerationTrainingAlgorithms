source("algos_Ai.R")
source("func_Ai.R")
source("plots.R")
theme_set(theme_bw())
library(MASS)
library(forcats)
load("iempaper.RData")
# save.image("iempaper.RData")

ni <- 6
n <- 1000
mu<-c(4,9)
mu0<-c(1,4)

da <- length(mu)
db <- 3

A <- matrix(sample.int(10, size = ni*ni, replace = TRUE), nrow = ni, ncol = da)
B <- matrix(sample.int(10, size = ni*ni, replace = TRUE), nrow = ni, ncol = db)
id <- rep(1:n,each=ni)
sigma <- 0.1*diag(ni)
omega <- 0.1*diag(db)

A <- vector("list", length=n)
B <- vector("list", length=n)
for (i in 1:n){
   Aj <- matrix(sample.int(10, size = ni*ni, replace = TRUE), nrow = ni, ncol = da)
   Bj <- matrix(sample.int(10, size = ni*ni, replace = TRUE), nrow = ni, ncol = db)
   A[[i]]<-Aj
   B[[i]]<-Bj
}

K<-3000
seed0=44444
nsim <- 3

col.names <- c("iteration", paste0("mu",1:da))
theta<-list(mu=mu)
theta0<-list(mu=mu0)

## EM
print('EM')
dem <- NULL
nbrem<-n
df.em <- vector("list", length=nsim)

nbriem1<-1
diem <- NULL
df.iem <- vector("list", length=nsim)

nbriem2<-n/2
diem50 <- NULL
df.iem50 <- vector("list", length=nsim)
for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  x <- NULL
  for (i in 1:n){
   xj <- mvrnorm(1, A[[i]]%*%mu,B[[i]]%*%omega%*%t(B[[i]])+sigma)
   x <- rbind(x,xj)
  }
  # print(head(x))
  df <- mixt.iem(x, theta0, K,A,B,sigma,omega,id,nbrem)
  ML <- df
  ML[1:(K+1),2:3]<- df[(K+1),2:3]
  df[,2:3] <- df[,2:3] - ML[,2:3]
  df$rep <- j
  dem <- rbind(dem,df)
  df$rep <- NULL
  df.em[[j]] <- df

  df <- mixt.iem(x, theta0, K,A,B,sigma,omega,id,nbriem1)
  ML <- df
  ML[1:(K+1),2:3]<- df[(K+1),2:3]
  df[,2:3] <- df[,2:3] - ML[,2:3]
  df$rep <- j
  diem <- rbind(diem,df)
  df$rep <- NULL
  df.iem[[j]] <- df

  df <- mixt.iem(x, theta0, K,A,B,sigma,omega,id,nbriem2)
  ML <- df
  ML[1:(K+1),2:3]<- df[(K+1),2:3]
  df[,2:3] <- df[,2:3] - ML[,2:3]
  df$rep <- j
  diem50 <- rbind(diem50,df)
  df$rep <- NULL
  df.iem50[[j]] <- df
}


# graphConvMC_new(dem, title="EM")
# graphConvMC_new(diem, title="IEM 1")
# graphConvMC_new(diem50, title="IEM 50")


dem[,2:3] <- dem[,2:3]^2
em <- NULL
em <- dem[dem$rep==1,]

if (nsim>2) {
   for (j in (2:nsim))
	{
	  em[,2:3] <- em[,2:3]+dem[dem$rep==j,2:3]
	}
}
em[,2:3] <- 1/nsim*em[,2:3]
em[,4]<-NULL


diem[,2:3] <- diem[,2:3]^2
iem1 <- NULL
iem1 <- diem[diem$rep==1,]

if (nsim>2) {
		for (j in (2:nsim))
	{
	  iem1[,2:3] <- iem1[,2:3]+diem[diem$rep==j,2:3]
	}
}

iem1[,2:3] <- 1/nsim*iem1[,2:3]
iem1[,4]<-NULL


diem50[,2:3] <- diem50[,2:3]^2
iem2 <- NULL
iem2 <- diem50[diem50$rep==1,]

if (nsim>2) {
    for (j in (2:nsim))
  {
    iem2[,2:3] <- iem2[,2:3]+diem50[diem50$rep==j,2:3]
  }
}

iem2[,2:3] <- 1/nsim*iem2[,2:3]
iem2[,4]<-NULL



em_scaled <- em
em_scaled$iteration = seq(0, n*K, by=n)
em_scaled <- em_scaled[rep(seq_len(nrow(em_scaled)), each=n),]

iem2_scaled <- iem2
iem2_scaled$iteration = seq(0, n/2*K, by=n/2)
iem2_scaled <- iem2_scaled[rep(seq_len(nrow(iem2_scaled)), each=n/2),]


iem1$algo <- 'IEM 1%'
iem2_scaled$algo <- 'IEM 50%'
em_scaled$algo <- 'EM'
variance <- NULL
variance <- rbind(em_scaled[0:K,],iem2_scaled[0:K,],iem1[0:K,])
colnames(variance) <- c("iteration","beta1","beta2","algo")
# graphConvMC2_new(variance, title="IEMs alpha=0.33",legend=TRUE)

error_1 <- function(df, title=NULL, ylim=NULL, legend=FALSE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df,aes(colour=df$algo ))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend ,size=1) +
      xlab("") + ylab(expression(paste(beta,"1"))) + theme_bw() +scale_linetype_manual(values = c("solid","longdash","dotted"))+scale_colour_manual(values = c('black','black','black')) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=10, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=10, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold",size=10)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}

error_2 <- function(df, title=NULL, ylim=NULL, legend=FALSE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df,aes(colour=df$algo ))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend ,size=1) +
      xlab("") +ylab(expression(paste(beta,"2"))) + theme_bw() +scale_linetype_manual(values = c("solid","longdash","dotted"))+scale_colour_manual(values = c('black','black','black')) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=10, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=10, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold",size=10)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}



graphConvMC_beta1_icml <- function(df,df2,df3, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  df2$algo <- as.factor(df2$algo)
  df3$algo <- as.factor(df3$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",linetype= "solid",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="black",linetype="longdash",size=1)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype="dotted",size=1)+
      xlab("") +ylab(expression(paste(beta,"1")))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=10, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=10, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=10)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}

graphConvMC_beta2_icml <- function(df,df2,df3, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  df2$algo <- as.factor(df2$algo)
  df3$algo <- as.factor(df3$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",linetype= "solid",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="black",linetype="longdash",size=1)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype="dotted",size=1)+
      xlab("") +ylab(expression(paste(beta,"2")))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=10, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=10, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=10)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}



variance <- NULL
variance <- rbind(em_scaled[0:4000,],iem2_scaled[0:4000,],iem1[0:4000,])
m <- error_1(variance[,c(1,2,4)])
n <- error_2(variance[,c(1,3,4)])
conv <- grid.arrange(m,n,ncol=2)

# plot_new3(variance,legend=FALSE)
m <- graphConvMC_beta1_icml(em_scaled[0:4000,c(1,2,4)],iem2_scaled[0:4000,c(1,2,4)],iem1[0:4000,c(1,2,4)])
n <- graphConvMC_beta2_icml(em_scaled[0:4000,c(1,3,4)],iem2_scaled[0:4000,c(1,3,4)],iem1[0:4000,c(1,3,4)])


############ FEW RUNS ################
K=40000
nsim=2
mus0<-list(c(1,5),c(3,7),c(2,6))
K2 <- 40000

# print('EM')
# dem <- NULL
# nbrem<-n
# df.em <- vector("list", length=nsim)

nbriem1<-1
diem <- NULL
df.iem <- vector("list", length=nsim)

# nbriem2<-n/2
# diem50 <- NULL
# df.iem50 <- vector("list", length=nsim)

for (j in (1:nsim))
{
  print(j)
  seed <- j*seed0
  set.seed(seed)
  theta0<-list(mu=mus0[[j]])
  x <- NULL
  for (i in 1:n){
   xj <- mvrnorm(1, A[[i]]%*%mu,B[[i]]%*%omega%*%t(B[[i]])+sigma)
   x <- rbind(x,xj)
  }
#   print('EM')
#   df <- mixt.iem(x, theta0, K,A,B,sigma,omega,id,nbrem)
#   df$rep <- j
#   df_scaled <- df
#   df_scaled$iteration = df_scaled$iteration*1
#   dem <- rbind(dem,df_scaled[1:(K2+1),])
#   df$rep <- NULL
#   df.em[[j]] <- df
# print('IEM 50')
#   df <- mixt.iem(x, theta0, K,A,B,sigma,omega,id,nbriem2)
#   df$rep <- j
#   df_scaled <- df
#   df_scaled$iteration = df_scaled$iteration*0.5
#   diem50 <- rbind(diem50,df_scaled[1:(K2+1),])
#   df$rep <- NULL
#   df.iem50[[j]] <- df
print('IEM 1/N')
  df <- mixt.iem(x, theta0, K,A,B,sigma,omega,id,nbriem1)
  df$rep <- j
  df_scaled <- df
  df_scaled$iteration = df_scaled$iteration*1/n
  diem <- rbind(diem,df_scaled[1:(K2+1),])
  df$rep <- NULL
  df.iem[[j]] <- df
}

dem$algo <- '100'
diem50$algo <- '50'
diem$algo <- '0.01'
colnames(dem) <- colnames(diem50)<- colnames(diem)<- c("iteration","beta1","beta2","rep","algo")



# run_1 <- function(df,df2,df3, title=NULL, ylim=NULL,legend=FALSE)
# {
#   G <- (ncol(df)-2)/3
#   df$rep <- as.factor(df$rep)
#   df2$rep <- as.factor(df2$rep)
#   df3$rep <- as.factor(df3$rep)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-3)
#   o <- c(0, 1, 2, 3, 4, 5, 6)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="red",linetype= "solid",size=0.5) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue",linetype="solid",size=0.5)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="green",linetype="solid",size=0.5)+
#     xlab("") + ylab(expression(paste(beta,"1")))+scale_x_continuous(limits = c(0, 70))+ theme_bw() +scale_linetype_manual(values = c("solid","solid","solid"))+scale_colour_manual(values = c('blue','green','red')) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
#                            size=10, angle=0),
#           axis.text.y = element_text(face="bold", color="black", 
#                            size=10, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black",face="bold", size=10)) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj
#   }
#   do.call("grid.arrange", c(graf, ncol=1, top=title))
# }



# run_2 <- function(df,df2,df3, title=NULL, ylim=NULL,legend=FALSE)
# {
#   G <- (ncol(df)-2)/3
#   df$rep <- as.factor(df$rep)
#   df2$rep <- as.factor(df2$rep)
#   df3$rep <- as.factor(df3$rep)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-3)
#   o <- c(0, 1, 2, 3, 4, 5, 6)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="red",linetype= "solid",size=0.5) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue",linetype="solid",size=0.5)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="green",linetype="solid",size=0.5)+
#     xlab("") + ylab(expression(paste(beta,"1")))+scale_x_continuous(limits = c(0, 70))+ theme_bw() +scale_linetype_manual(values = c("solid","solid","solid"))+scale_colour_manual(values = c('blue','green','red')) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
#                            size=10, angle=0),
#           axis.text.y = element_text(face="bold", color="black", 
#                            size=10, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black",face="bold", size=10)) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj
#   }
#   do.call("grid.arrange", c(graf, ncol=1, top=title))
# }


# # dem$algo <- 'EM'
# # diem50$algo <- 'MBEM50'
# # diem$algo <- 'IEM'
# a <- run_1(dem[,c(1,2,4)],diem50[,c(1,2,4)],diem[,c(1,2,4)])
# b <- run_2(dem[,c(1,3,4)],diem50[,c(1,3,4)],diem[,c(1,3,4)])
# grid.arrange(a,b,ncol=2)



seplot <- function(df,colname,title.labs, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  df$rep <- as.factor(df$rep)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  graf <- ggplot(df)+geom_line(aes(iteration,value,by=value,colour = df$algo,linetype=df$rep),show.legend = legend,size=1)+guides(linetype=FALSE,size=FALSE)+labs(colour=title.labs) +
  xlab("passes")+ ylab(colname)+scale_linetype_manual(values=c("solid", "solid"))  + theme_bw() +scale_x_continuous(limits = c(0, 70))+ theme(legend.position = c(0.8, 0.4)) + guides(color = guide_legend(override.aes = list(size=5)))
   theme(legend.text=element_text(size=20),legend.title=element_text(size=20))+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(color="black", 
                           size=20, angle=0),
          axis.text.y = element_text(color="black", 
                           size=20, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", size=20)) + theme(aspect.ratio=1)
  grid.arrange(graf)
  # do.call("grid.arrange", c(graf, ncol=1, top=title))
}


seplot <- function(df,colname,title.labs, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  df$rep <- as.factor(df$rep)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  graf <- ggplot(df)+geom_line(aes(iteration,value,by=value,colour = df$algo,linetype=df$rep),show.legend = legend,size=1)+guides(linetype=FALSE,size=FALSE)+labs(colour=title.labs) +
  xlab("passes")+ ylab(colname)+scale_linetype_manual(values=c("solid", "solid"))  +scale_x_continuous(limits = c(0, 70))  + theme_bw() + theme(
        
        panel.background = element_rect(colour = "grey", size=1),legend.position = c(0.7, 0.6)) + guides(color = guide_legend(override.aes = list(size=9)))+
   theme(legend.text=element_text(size=24),legend.title=element_text(size=24))+ theme(panel.border = element_blank() ,axis.text.x = element_text(color="black", 
                           size=20, angle=0),
          axis.text.y = element_text(color="black", 
                           size=20, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", size=20)) + theme(aspect.ratio=1)
  grid.arrange(graf)
  # do.call("grid.arrange", c(graf, ncol=1, top=title))
}


dem$algo <- 'p=N'
diem50$algo <- 'p=N/2'
diem$algo <- 'p=1'
colnames(dem) <- colnames(diem50)<- colnames(diem)<- c("iteration","beta1","beta2","rep","algo")


comparison <- 0
comparison <- rbind(dem[,c(1,2,4,5)],diem50[,c(1,2,4,5)],diem[,c(1,2,4,5)])
comparison$algo <- factor(comparison$algo, levels=c('p=1','p=N/2','p=N'), ordered = TRUE)
var <- melt(comparison, id.var = c('iteration','algo','rep'), na.rm = TRUE)

beta0 <- seplot(var,expression(paste(theta,"1")),title.labs = 'batch size',title="comparison",legend=FALSE)

comparison <- 0
comparison <- rbind(dem[,c(1,3,4,5)],diem50[,c(1,3,4,5)],diem[,c(1,3,4,5)])
comparison$algo <- factor(comparison$algo, levels=c('p=1','p=N/2','p=N'), ordered = TRUE)
var <- melt(comparison, id.var = c('iteration','algo','rep'), na.rm = TRUE)

gamma <- seplot(var,expression(paste(theta,"2")),title.labs = 'batch size', title="comparison",legend=TRUE)

final <- grid.arrange(beta0,gamma,ncol=2)

# ggsave(final,file="iems.pdf", width = 900, height = 450, units = "mm")