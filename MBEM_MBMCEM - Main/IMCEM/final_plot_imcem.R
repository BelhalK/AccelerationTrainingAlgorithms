load("final_cat.RData")
# save.image("final_cat2.RData")

setwd("/Users/karimimohammedbelhal/Desktop/imcem_logistic")
source('plots.R') 

require(ggplot2)
require(gridExtra)
require(reshape2)






comparison <- 0

comparison <- rbind(
  theo_mix25_scaledbis[1:end,c(1,2,8)],
  theo_mix50_scaledbis_longer[,c(1,2,8)],
  theo_mix35_scaledbis[1:end,c(1,2,8)],theo_mix85_scaledbis[1:end,c(1,2,8)],theo_ref_scaledbis[1:end,c(1,2,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

beta0 <- seplot(var,expression(paste(beta,"1")), title="comparison",legend=FALSE)

comparison <- 0
comparison <- rbind(
  theo_mix25_scaledbis[1:end,c(1,3,8)],
  theo_mix50_scaledbis_longer[,c(1,3,8)],
  theo_mix35_scaledbis[1:end,c(1,3,8)],theo_mix85_scaledbis[1:end,c(1,3,8)],theo_ref_scaledbis[1:end,c(1,3,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

gamma <- seplot(var,expression(paste(beta,"2")), title="comparison",legend=FALSE)

comparison <- 0
comparison <- rbind(
  theo_mix25_scaledbis[1:end,c(1,4,8)],
  theo_mix50_scaledbis_longer[,c(1,4,8)],
  theo_mix35_scaledbis[1:end,c(1,4,8)],theo_mix85_scaledbis[1:end,c(1,4,8)],theo_ref_scaledbis[1:end,c(1,4,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

delta <- seplot(var,expression(paste(beta,"3")), title="comparison",legend=TRUE)


final <- grid.arrange(beta0,gamma,delta,ncol=3)




seplotgamma <- function(df,colname, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  scaleFUN <- function(x) sprintf("%.1f", x)
  graf <- vector("list", ncol(df)-2)
  graf <- ggplot(df)+geom_line(aes(iterations,value,by=value,colour = df$algo),show.legend = legend)+labs(colour='batch size') +
  xlab("passes")+ ylab(colname)+  scale_y_continuous(breaks = c(-0.5,-0.4,-0.3))  + theme_bw() + theme(
        
        panel.background = element_rect(colour = "grey", size=1),legend.position = c(0.8, 0.6)) + guides(color = guide_legend(override.aes = list(size=5)))+
   theme(legend.text=element_text(size=20),legend.title=element_text(size=20))+ theme(panel.border = element_blank() ,axis.text.x = element_text(color="black", 
                           size=20, angle=0),
          axis.text.y = element_text(color="black", 
                           size=20, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", size=20)) + theme(aspect.ratio=1)
  grid.arrange(graf)
  # do.call("grid.arrange", c(graf, ncol=1, top=title))
}

seplot <- function(df,colname, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  scaleFUN <- function(x) sprintf("%.1f", x)
  graf <- vector("list", ncol(df)-2)
  graf <- ggplot(df)+geom_line(aes(iterations,value,by=value,colour = df$algo),show.legend = legend)+labs(colour='batch size') +
  xlab("passes")+ ylab(colname) + theme_bw() + theme(
        
        panel.background = element_rect(colour = "grey", size=1),legend.position = c(0.8, 0.6)) + guides(color = guide_legend(override.aes = list(size=5)))+
   theme(legend.text=element_text(size=20),legend.title=element_text(size=20))+ theme(panel.border = element_blank() ,axis.text.x = element_text(color="black", 
                           size=20, angle=0),
          axis.text.y = element_text(color="black", 
                           size=20, angle=0))+theme(axis.title = element_text(family = "Trebuchet MS", color="black", size=20)) + theme(aspect.ratio=1)
  grid.arrange(graf)
  # do.call("grid.arrange", c(graf, ncol=1, top=title))
}

#fixed effects
comparison <- 0

comparison <- rbind(
  theo_mix25_scaledbis[10:end,c(1,2,8)],
  theo_mix50_scaledbis_longer[10:nrow(theo_mix50_scaledbis_longer),c(1,2,8)],
  theo_mix35_scaledbis[10:end,c(1,2,8)],theo_mix85_scaledbis[10:end,c(1,2,8)],theo_ref_scaledbis[10:end,c(1,2,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

beta0 <- seplot(var,expression(paste(beta,"1")), title="comparison",legend=FALSE)

comparison <- 0
comparison <- rbind(
  theo_mix25_scaledbis[10:end,c(1,3,8)],
  theo_mix50_scaledbis_longer[10:nrow(theo_mix50_scaledbis_longer),c(1,3,8)],
  theo_mix35_scaledbis[10:end,c(1,3,8)],theo_mix85_scaledbis[10:end,c(1,3,8)],theo_ref_scaledbis[10:end,c(1,3,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

gamma <- seplotgamma(var,expression(paste(beta,"2")), title="comparison",legend=FALSE)

comparison <- 0
comparison <- rbind(
  theo_mix25_scaledbis[10:end,c(1,4,8)],
  theo_mix50_scaledbis_longer[10:nrow(theo_mix50_scaledbis_longer),c(1,4,8)],
  theo_mix35_scaledbis[10:end,c(1,4,8)],theo_mix85_scaledbis[10:end,c(1,4,8)],theo_ref_scaledbis[10:end,c(1,4,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

delta <- seplot(var,expression(paste(beta,"3")), title="comparison",legend=TRUE)


final <- grid.arrange(beta0,gamma,delta,ncol=3)



#variances
comparison <- 0

comparison <- rbind(
  theo_mix25_scaledbis[10:end,c(1,5,8)],
  theo_mix50_scaledbis_longer[10:nrow(theo_mix50_scaledbis_longer),c(1,5,8)],
  theo_mix35_scaledbis[10:end,c(1,5,8)],theo_mix85_scaledbis[10:end,c(1,5,8)],theo_ref_scaledbis[10:end,c(1,5,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

beta0 <- seplot(var,expression(paste(beta,"1")), title="comparison",legend=FALSE)

comparison <- 0
comparison <- rbind(
  theo_mix25_scaledbis[10:end,c(1,6,8)],
  theo_mix50_scaledbis_longer[10:nrow(theo_mix50_scaledbis_longer),c(1,6,8)],
  theo_mix35_scaledbis[10:end,c(1,6,8)],theo_mix85_scaledbis[10:end,c(1,6,8)],theo_ref_scaledbis[10:end,c(1,6,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

gamma <- seplotgamma(var,expression(paste(beta,"2")), title="comparison",legend=FALSE)

comparison <- 0
comparison <- rbind(
  theo_mix25_scaledbis[10:end,c(1,7,8)],
  theo_mix50_scaledbis_longer[10:nrow(theo_mix50_scaledbis_longer),c(1,7,8)],
  theo_mix35_scaledbis[10:end,c(1,7,8)],theo_mix85_scaledbis[10:end,c(1,7,8)],theo_ref_scaledbis[10:end,c(1,7,8)])

var <- melt(comparison, id.var = c('iterations','algo'), na.rm = TRUE)

delta <- seplot(var,expression(paste(beta,"3")), title="comparison",legend=TRUE)


final <- grid.arrange(beta0,gamma,delta,ncol=3)


graphConvMC_5(theo_ref,theo_mix25bis,theo_mix,theo_mix85,theo_mix35)
graphConvMC_5(theo_ref,theo_mix25,theo_mix,theo_mix85,theo_mix25)
