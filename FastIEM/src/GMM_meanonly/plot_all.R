require(ggplot2)
require(gridExtra)
require(reshape2)
library(rlist)

source("utils/algos.R")
source("utils/func.R")
source("utils/plots.R")
theme_set(theme_bw())
options(digits = 22)


datasizes <- seq(1000, 106000, 5000)


load("RData_sep/emiter.RData")
load("RData_sep/iemiter.RData")
load("RData_sep/iemseqiter.RData")
load("RData_sep/oemiter.RData")
load("RData_sep/oemvriter.RData")
load("RData_sep/sagaiter.RData")





length(datasizes)



# eqiem = function(x){x+100}
# eqsaga = function(x){x**(2/3)}


# for (precision in c(1e-3,1e-4,1e-5)){
#   emindex  <- iemseqindex <- iemindex  <- oemvrindex <- sagaindex  <- c()

#   for (i in (1:length(datasizes))){
#   	emindex[i] <- datasizes[[i]]*which(emiter[[i]][,c(4)] < precision)[1]
#     iemindex[i] <- which(iemiter[[i]][,c(4)] < precision)[1]
#     iemseqindex[i] <- which(iemseqiter[[i]][,c(4)] < precision)[1]
#     sagaindex[i] <- which(sagaiter[[i]][,c(4)] < precision)[1]
#     oemvrindex[i] <- which(oemvriter[[i]][,c(4)] < precision)[1]
#   }

#   x  <- datasizes
#   y1 <- iemindex
#   y2 <- iemseqindex
#   y3 <- sagaindex
#   y4 <- oemvrindex
#   y5 <- emindex
#   ytwothird <- eqsaga(datasizes)
#   ylin <- eqiem(datasizes)
#   df <- data.frame(x,y1,y2,y4,y5,ytwothird,ylin)

#   print(ggplot(df, aes(x),show.legend = TRUE) +                    
#     geom_line(aes(y=y1), colour="red") +
#     geom_line(aes(y=y3), colour="pink") +
#   geom_line(aes(y=y2), colour="blue") +
#   geom_line(aes(y=y4), colour="purple") +
#   geom_line(aes(y=y5), colour="green") +
#   geom_line(aes(y=ytwothird), colour="black", linetype= "dotted") +
#   geom_line(aes(y=ylin), colour="brown", linetype= "dotted") +
#   scale_y_log10()  + scale_x_log10() +
#   xlab("Dataset size") + ylab("Iteration")  +
#   ggtitle(precision))

# }







eqiem = function(x){x*5}
eqsaga = function(x){(x*10)**(2/3)}

precision = 1e-3
emindex  <- iemseqindex <- iemindex  <- oemvrindex <- sagaindex  <- c()

for (i in (1:length(datasizes))){
  emindex[i] <- datasizes[[i]]*which(emiter[[i]][,c(4)] < precision)[1]
  iemindex[i] <- which(iemiter[[i]][,c(4)] < precision)[1]
  iemseqindex[i] <- which(iemseqiter[[i]][,c(4)] < precision)[1]
  sagaindex[i] <- which(sagaiter[[i]][,c(4)] < precision)[1]
  oemvrindex[i] <- which(oemvriter[[i]][,c(4)] < precision)[1]
}

x  <- datasizes
y1 <- data.frame(iteration=datasizes, prec = iemindex, algo = "IEM")
y2 <- data.frame(iteration=datasizes, prec = iemseqindex, algo = "IEM")
y3 <- data.frame(iteration=datasizes, prec = sagaindex, algo = "FI-EM")
y4 <- data.frame(iteration=datasizes, prec = oemvrindex, algo = "SVR-EM")
y5 <- data.frame(iteration=datasizes, prec = emindex, algo = "EM")

ytwothird <- data.frame(iteration=datasizes, prec = eqsaga(datasizes), algo = "f(n) = n^(2/3)")
ylin <- data.frame(iteration=datasizes, prec = eqiem(datasizes), algo = "f(n) = n")

sublin = eqsaga(datasizes)
lin = eqiem(datasizes)
list(sublin)[[1]]
write.table(list(sublin)[[1]], 'notebooks/test.txt')

df <- rbind(y2,y3,y4,y5,ytwothird,ylin)
colnames(df) <- c("iteration","prec","algo")

write.csv(df, file = "notebooks/rates.csv")
# df.m <- melt(df, id.var = c("iteration","algo"))

# testdf <- df.m[,c(1,4,2)]


plotn <- function(df, title=NULL, ylim=NULL, legend=TRUE)
{
  G <- (ncol(df)-2)/3
  df$algo <- as.factor(df$algo)
  ylim <-rep(ylim,each=2)
  graf <- vector("list", ncol(df)-2)
  o <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  for (j in (2:(ncol(df)-1)))
  {
    grafj <- ggplot(df,aes(colour=df$algo ,linetype=df$algo))+geom_line(aes_string(df[,1],df[,j],by=df[,ncol(df)]),show.legend = legend) + 
    scale_linetype_manual(values = c("solid","solid","solid","solid","dotted","dotted")) +
      xlab("epochs")+ scale_y_log10()  + scale_x_log10() +xlab("Problem size n") + ylab("")+
      guides(color = guide_legend(override.aes = list(size = 2))) +
      theme(axis.text.x = element_text(face="bold", color="black", 
                           size=20, angle=0),
          axis.text.y = element_text(face="bold", color="black", 
                           size=20, angle=0),axis.title = element_text( color="black", face="bold",size=20),legend.text=element_text(size=20),
          legend.title = element_blank(),legend.position = c(0.2, 0.8))
    if (!is.null(ylim))
      grafj <- grafj + ylim(ylim[j-1]*c(-1,1))
    graf[[o[j]]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol=1))
}

plotn(df)


save <- plotn(df)
ggsave(save,file="/Users/karimimohammedbelhal/Desktop/rates.pdf", width = 200, height = 200, units = "mm")




# em_ep <- emiter
# for (i in (1:length(datasizes))){
#   em_ep[[i]]$iteration <- datasizes[i]*emiter[[i]]$iteration
# }



# i = 15
# start = 0
# end = 20000
# variance <- rbind(emiter[[i]][start:end,c(1,4,8)],
#                   iemseqiter[[i]][start:end,c(1,4,8)],
#                   sagaiter[[i]][start:end,c(1,4,8)],
#                    oemvriter[[i]][start:end,c(1,4,8)])

# graphConvMC2_new(variance, title="IEMs GMM 1e5",legend=TRUE)


# length(oemvrl)
# curveoemvr <- oemvriter
# curvesaga <- sagaiter
# curveem <- emiter
# curveiem <- iemiter
# curveiemseq <- iemseqiter

# for (i in (1:length(datasizes))){
#   curveoemvr[[i]]$size <- datasizes[i]
#   curvesaga[[i]]$size <- datasizes[i]
#   curveem[[i]]$size <- datasizes[i]
#   curveiem[[i]]$size <- datasizes[i]
#   curveiemseq[[i]]$size <- datasizes[i]
# }



# epochs
# start =1
# end = 100000
# plot.oemvr <- NULL
# plot.saga <- NULL
# for (i in (1:length(datasizes))){
#   plot.oemvr <- rbind(plot.oemvr, curveoemvr[[i]][start:end,c(1,4,9)])
#   plot.saga <- rbind(plot.saga, curvesaga[[i]][start:end,c(1,4,9)])

# }

# plotagainstn(plot.oemvr, title="SAGA-EM for different dataset size",legend=TRUE)
# plotagainstn(plot.saga, title="SAGA-EM for different dataset size",legend=TRUE)




# ### EPOCH WISE

# ### PER EPOCH
# eml <- ieml <- iemseql <- oeml <- oemvrl <- sagal <- list()


# em_ep <- emiter
# for (i in (1:length(datasizes))){
#   em_ep[[i]]$iteration <- datasizes[i]*emiter[[i]]$iteration
# }


# for (i in (1:length(datasizes))){
#   n <- datasizes[i]
#   K = 20*n
#   epochs = seq(1, K, by=n)
#   iem_ep <- iemiter[[i]][epochs,]
#   iem_ep$iteration <- 1:(K/n)
#   iemseq_ep <- iemseqiter[[i]][epochs,]
#   iemseq_ep$iteration <- 1:(K/n)
#   oem_ep <- oemiter[[i]][epochs,]
#   oem_ep$iteration <- 1:(K/n)
#   oemvr_ep <- oemvriter[[i]][epochs,]
#   oemvr_ep$iteration <- 1:(K/n)
#   saga_ep <- sagaiter[[i]][epochs,]
#   saga_ep$iteration <- 1:(K/n)


#   ieml[[i]] <- iem_ep
#   iemseql[[i]] <- iemseq_ep
#   oeml[[i]] <- oem_ep
#   oemvrl[[i]] <- oemvr_ep
#   sagal[[i]] <- saga_ep
# }


# curveoemvr <- oemvrl
# curvesaga <- sagal
# curveiem <- ieml
# curveiemseq <- iemseql

# for (i in (1:length(datasizes))){
#   curveoemvr[[i]]$size <- datasizes[i]
#   curvesaga[[i]]$size <- datasizes[i]
#   curveiem[[i]]$size <- datasizes[i]
#   curveiemseq[[i]]$size <- datasizes[i]
# }



# epochs
# start =1
# end = 10
# plot.iemseq <- NULL
# plot.oem <- NULL
# plot.oemvr <- NULL
# plot.saga <- NULL
# for (i in (1:length(datasizes))){
#   plot.iemseq <- rbind(plot.iemseq, curveiemseq[[i]][start:end,c(1,4,9)])
#   plot.oemvr <- rbind(plot.oemvr, curveoemvr[[i]][start:end,c(1,4,9)])
#   plot.saga <- rbind(plot.saga, curvesaga[[i]][start:end,c(1,4,9)])

# }

# plotagainstn(plot.oemvr, title="OEM-vr for different dataset size",legend=TRUE)
# plotagainstn(plot.saga, title="SAGA-EM for different dataset size",legend=TRUE)
# plotagainstn(plot.iemseq, title="IEM Seq for different dataset size",legend=TRUE)




# ### JUST ONE RUN



# start =1
# end = 100000

# i = 15
# variance <- rbind(oemvriter[[i]][start:end,c(1,4,8)],
#                   iemiter[[i]][start:end,c(1,4,8)],
#                   oemiter[[i]][start:end,c(1,4,8)],
#                   iemseqiter[[i]][start:end,c(1,4,8)],
#                   sagaiter[[i]][start:end,c(1,4,8)])

# graphConvMC2_new(variance, title="IEMs GMM 1e5",legend=TRUE)


