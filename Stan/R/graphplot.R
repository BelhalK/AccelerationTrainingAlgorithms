
graphConvMC_twokernels <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df))))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+
      xlab("iteration")+ ylab(names(df[j])) + theme_bw()
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}


plotmcmc <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df))))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue",linetype="solid")+
      xlab("iteration")+ ylab(names(df[j])) + theme_bw()
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=ncol(df)-1, top=title))
}



plotquantile <- function(df,df2, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=1)+
      xlab("")+scale_x_log10()+ theme_bw() +ylab(names(df[j]))+ theme(axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=15, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=20)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}


plotquantile3 <- function(df,df2,df3, title=NULL, ylim=NULL)
{
 G <- (ncol(df)-2)/3
  df$quantile <- as.factor(df$quantile)
  df2$quantile <- as.factor(df2$quantile)
  df3$quantile <- as.factor(df3$quantile)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=1) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=1)+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype = 2,size=1)+
      xlab("")+scale_x_log10()+ theme_bw() +ylab(names(df[j]))+ theme(axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=15, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=20)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}

plotconv <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df))))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+
      xlab("iteration")+scale_x_log10()+ ylab(names(df[j])) + theme_bw()
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=2, top=title))
}




plotconv3 <- function(df,df2,df3, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df))))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red")+
      xlab("iteration")+ ylab(names(df[j])) + theme_bw()
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=2, top=title))
}

plotconv4 <- function(df,df2,df3,df4, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df))))
  {
    grafj <- ggplot(df)+geom_point(aes_string(df[,1],df[,j],by=df[,ncol(df)]),shape=1)+
    geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]))+
    geom_point(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+
    geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+
    geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="black",linetype="dotted")+
    geom_line(aes_string(df4[,1],df4[,j],by=df4[,ncol(df4)]),colour="blue",linetype="dotted")+
      xlab("iteration")+scale_x_log10()+ ylab(names(df[j])) + theme_bw()
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=2, top=title))
}


graphConvMC_new <- function(df, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +
      xlab("iteration") + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=ncol(df)-2, top=title))
}

graphConvMC_diff <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red")+
      xlab("iteration") +scale_x_log10()+ ylab(names(df[j]))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=2, top=title))
}



graphConvMC_diffzoom <- function(df,df2,df3, K1,end,title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red")+
      xlab("iteration")+scale_x_continuous(limits = c(K1, end))+scale_y_continuous(limits = c(df[end,j]-0.2, df[end,j]+0.2))+ ylab(names(df[j]))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=2, top=title))
}