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
  do.call("grid.arrange", c(graf, ncol=ncol(df)-2, top=title))
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
      xlab("iteration") + ylab(names(df[j])) + theme_bw()
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}


graphConvMC_threekernels <- function(df,df2,df3, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue")+geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red")+
      xlab("iteration") + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=3, top=title))
}


graphConvMC3_new <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  graf <- ggplot(df,aes(colour=df$algo ))+geom_line(aes(iteration,value,by=value),show.legend = legend) +
  xlab("iteration") + ylab('value') + facet_wrap(~variable,scales = "free_y") #+ coord_trans(x = "log10")
  grid.arrange(graf)
  # do.call("grid.arrange", c(graf, ncol=1, top=title))
}




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