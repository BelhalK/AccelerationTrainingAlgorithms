

plot2 <- function(df,df2, title=NULL, ylim=NULL)
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

plot3 <- function(df,df2,df3, title=NULL, ylim=NULL)
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
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}


plot5 <- function(df,df2,df3,df4,df5, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df))))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +
    geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+
    geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red")+
    geom_line(aes_string(df4[,1],df4[,j],by=df4[,ncol(df4)]),colour="pink")+
    geom_line(aes_string(df5[,1],df5[,j],by=df5[,ncol(df5)]),colour="yellow")+
      xlab("iteration")+ ylab(names(df[j])) + theme_bw()
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}



diff <- function(df,df2,df3, title=NULL, ylim=NULL)
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
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",size=1) +
    geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue",linetype = 2,size=1)+
    geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red",linetype = 2,size=1)+
      xlab("") +scale_x_log10(breaks= c(10,100,200))+ ylab(names(df[j])) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text( color="black", 
                           size=10, angle=0),
          axis.text.y = element_text( color="black", 
                           size=10, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black",  size=10)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}