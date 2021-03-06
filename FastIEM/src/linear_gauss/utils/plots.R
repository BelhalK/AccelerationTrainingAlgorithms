require(ggplot2)
require(gridExtra)
require(reshape2)



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



# plot_new <- function(df, title=NULL, ylim=NULL, legend=TRUE)
# {
#   G <- (ncol(df)-2)/3
#   df$algo <- as.factor(df$algo)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df,aes(colour=df$algo ))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend ,size=1) +
#       xlab("") + ylab(names(df[j])) + theme_bw() +scale_linetype_manual(values = c("dashed","solid"))+scale_colour_manual(values = c('blue','red')) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
#                            size=10, angle=0),
#           axis.text.y = element_text(face="bold", color="black", 
#                            size=10, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", size=13)) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj
#   }
#   do.call("grid.arrange", c(graf, ncol=1, top=title))
# }

# plot_new3 <- function(df, title=NULL, ylim=NULL, legend=TRUE)
# {
#   G <- (ncol(df)-2)/3
#   df$algo <- as.factor(df$algo)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df,aes(colour=df$algo ))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend ,size=1) +
#       xlab("") + ylab(names(df[j])) + theme_bw() +scale_linetype_manual(values = c("solid","solid","solid"))+scale_colour_manual(values = c('blue','green','red')) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
#                            size=15, angle=0),
#           axis.text.y = element_text(face="bold", color="black", 
#                            size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold",size=20)) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj
#   }
#   do.call("grid.arrange", c(graf, ncol=2, top=title))
# }

# graphConvMC_diff4 <- function(df,df2, title=NULL, ylim=NULL)
# {
#   G <- (ncol(df)-2)/3
#   df$individual <- as.factor(df$individual)
#   df2$individual <- as.factor(df2$individual)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue",size=2) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype = 2,size=2)+
#       xlab("") +scale_x_log10()+ ylab(names(df[j]))  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
#                            size=20, angle=0),
#           axis.text.y = element_text(face="bold", color="black", 
#                            size=20, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=30)) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj

#   }
#   do.call("grid.arrange", c(graf, ncol=1, top=title))
# }


# graphConvMC_new <- function(df, title=NULL, ylim=NULL)
# {
#   G <- (ncol(df)-2)/3
#   df$rep <- as.factor(df$rep)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +
#       xlab("iteration") + ylab(names(df[j])) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj
#   }
#   do.call("grid.arrange", c(graf, ncol=1, top=title))
# }

# graphConvMCdf2_new <- function(df,df2, title=NULL, ylim=NULL)
# {
#   G <- (ncol(df)-2)/3
#   df$rep <- as.factor(df$rep)
#   df2$rep <- as.factor(df2$rep)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)])) +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="blue",size=1) +
#       xlab("iteration") + ylab(names(df[j])) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj
#   }
#   do.call("grid.arrange", c(graf, ncol=1, top=title))
# }

# graphConvMCdf3_new <- function(df,df2,df3, title=NULL, ylim=NULL)
# {
#   G <- (ncol(df)-2)/3
#   df$rep <- as.factor(df$rep)
#   df2$rep <- as.factor(df2$rep)
#   df3$rep <- as.factor(df3$rep)
#   df$algo <- as.factor(df$algo)
#   df2$algo <- as.factor(df2$algo)
#   df3$algo <- as.factor(df3$algo)
#   ylim <-rep(ylim,each=2)
#   graf <- vector("list", ncol(df)-2)
#   o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#   for (j in (2:(ncol(df)-1)))
#   {
#     grafj <- ggplot(df)+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),colour="blue") +geom_line(aes_string(df2[,1],df2[,j],by=df2[,ncol(df2)]),colour="red",linetype=1) +geom_line(aes_string(df3[,1],df3[,j],by=df3[,ncol(df3)]),colour="green",linetype=1)+
#       xlab("") + ylab(names(df[j])) + theme_bw() +scale_linetype_manual(values = c("solid","solid","solid"))+scale_colour_manual(values = c('blue','green','red')) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="black", 
#                            size=15, angle=0),
#           axis.text.y = element_text(face="bold", color="black", 
#                            size=15, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black",face="bold", size=20)) 
#     if (!is.null(ylim))
#       grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
#     graf[[o[j]]] <- grafj
#   }
#   do.call("grid.arrange", c(graf, ncol=2, top=title))
# }


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
      xlab("epoch") +scale_y_log10() + ylab(names(df[j])) 
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=1, top=title))
}
