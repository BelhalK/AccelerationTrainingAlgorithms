require(ggplot2)
require(gridExtra)
require(reshape2)



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
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}



graphConvMC_twokernels <- function(df,df2, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+
      xlab("iteration") + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}

