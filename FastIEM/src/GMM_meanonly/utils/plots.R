

mixt.ident <- function(df)
{
  G <- (ncol(df)-1)/3
  K <- nrow(df)
  mu.final <- as.numeric(as.character(df[K,(G+2):(2*G+1)]))
  ind <- sort.int(mu.final, index.return=TRUE)$ix
  df[,2:(G+1)] <- df[,(G-1+ind)]
  df[,(G+2):(2*G+1)] <- df[,(2*G-1+ind)]
  df[,(2*G+2):(3*G+1)] <- df[,(3*G-1+ind)]
  return(df)
}

mixt.ident3 <- function(df)
{
  G <- (ncol(df)-1)/3
  K <- nrow(df)
  mu.final <- as.numeric(as.character(df[K,(G+2):(2*G+1)]))
  ind <- sort.int(mu.final, index.return=TRUE)$ix
  df[,2:(G+1)] <- df[,(G-2+ind)]
  df[,(G+2):(2*G+1)] <- df[,(2*G-2+ind)]
  df[,(2*G+2):(3*G+1)] <- df[,(3*G-2+ind)]
  return(df)
}



graphConvMC_new <- function(df, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$rep <- as.factor(df$rep)
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


graphConvMC2_new <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df,aes(colour=df$algo ))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend) +
      xlab("epochs")+ scale_y_log10() +ylab("")+scale_x_continuous(breaks = round(seq(1, 14, by = 2),1))   + 
      guides(color = guide_legend(override.aes = list(size = 2))) +
      theme(axis.text.x = element_text(face="bold", color="black", 
                           size=20, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=20, angle=0),axis.title = element_text( color="black", face="bold",size=20),legend.text=element_text(size=20),
          legend.title = element_blank(),legend.position = c(0.2, 0.2))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=1))
}



graphConvMC <- function(df,df2,df3,df4, title=NULL, ylim=NULL)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  df2$algo <- as.factor(df2$algo)
  df3$algo <- as.factor(df3$algo)
  df4$algo <- as.factor(df4$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="black",linetype= "solid",size=2)+
    geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="black",linetype="longdash",size=2)+
    geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="red",linetype="dotted",size=2)+
    geom_line(aes_string(df4[,1],df4[,j],by=df4[,ncol(df4)]),colour="blue",linetype="dotted",size=2)+
      xlab("") +ylab(expression(paste(beta,"1")))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
                           size=30, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=30, angle=0))+theme(axis.title = element_text(color="black", face="bold", size=30)) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj

  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}




plotagainstn <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$size <- as.factor(df$size)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df,aes(colour=df$size ))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend) +
      xlab("epochs")+ scale_y_log10()  + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}
