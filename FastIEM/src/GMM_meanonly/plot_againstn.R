require(ggplot2)
require(gridExtra)
require(reshape2)
library(rlist)

source("utils/algos.R")
source("utils/func.R")
source("utils/plots.R")
theme_set(theme_bw())
options(digits = 22)

# save.image("RData/precisionagainstn.RData")
# load("RData/precisionagainstn_VM.RData")
load("RData_VM/precisionagainstn_VM.RData")
# load("RData_VM/precisionagainstn_VM_onlyoemvriem.RData")
load("RData/iem.RData")
load("RData/iemseq.RData")
load("RData/oemvr.RData")

load("RData/ieml.RData")
load("RData/iemseql.RData")
load("RData/oemvrl.RData")



# iemiternew <- iemiter
# iemseqiternew <- iemseqiter
# oemvriternew <- oemvriter


length(datasizes)


em_ep <- eml
for (i in (1:length(datasizes))){
	em_ep[[i]]$iteration <- datasizes[i]*eml[[i]]$iteration
}

eqiem = function(x){x+100}
eqsaga = function(x){x**(2/3)}

# x  <- datasizes
# y1 <- eqiem(datasizes)
# y2 <- eqsaga(datasizes)
# df <- data.frame(x,y1,y2)

# ggplot(df, aes(x),show.legend = TRUE) +                    
#   geom_line(aes(y=y1), colour="red") +  
#   geom_line(aes(y=y2), colour="green") +
#   xlab("Dataset size") + ylab("Epoch")  +
#   ggtitle(precision)


# plot(eqsaga(1000:100000), type='l')



for (precision in c(1e-2,1e-3,1e-4,1e-5)){
  emindex  <- iemseqindex <- iemindex  <- oemvrindex <- sagaindex  <- c()

  for (i in (1:length(datasizes))){
  	emindex[i] <- datasizes[[i]]*which(emiter[[i]][,c(4)] < precision)[1]
    iemindex[i] <- which(iemiternew[[i]][,c(4)] < precision)[1]
    iemseqindex[i] <- which(iemseqiternew[[i]][,c(4)] < precision)[1]
    sagaindex[i] <- which(sagaiter[[i]][,c(4)] < precision)[1]
    oemvrindex[i] <- which(oemvriternew[[i]][,c(4)] < precision)[1]
  }

  x  <- datasizes
  y1 <- iemindex
  y2 <- iemseqindex
  y3 <- sagaindex
  y4 <- oemvrindex
  y5 <- emindex
  ytwothird <- eqsaga(datasizes)
  ylin <- eqiem(datasizes)
  df <- data.frame(x,y1,y2,y4,y5,ytwothird,ylin)

  print(ggplot(df, aes(x),show.legend = TRUE) +                    
    geom_line(aes(y=y1), colour="red") +
    geom_line(aes(y=y3), colour="pink") +
  geom_line(aes(y=y2), colour="blue") +
  geom_line(aes(y=y4), colour="purple") +
  geom_line(aes(y=y5), colour="green") +
  geom_line(aes(y=ytwothird), colour="black", linetype= "dotted") +
  geom_line(aes(y=ylin), colour="brown", linetype= "dotted") +
  scale_y_log10()  + scale_x_log10() +
  xlab("Dataset size") + ylab("Iteration")  +
  ggtitle(precision))

}




for (precision in c(1e-2,1e-3,1e-4,1e-5)){
  iemseqindex <- iemindex  <- oemvrindex <- c()

  for (i in (1:length(datasizes))){
    iemindex[i] <- which(iemiter[[i]][,c(4)] < precision)[1]
    iemseqindex[i] <- which(iemseqiter[[i]][,c(4)] < precision)[1]
    # sagaindex[i] <- which(sagaiter[[i]][,c(4)] < precision)[1]
    oemvrindex[i] <- which(oemvriter[[i]][,c(4)] < precision)[1]
  }

  x  <- datasizes
  y1 <- iemindex
  y2 <- iemseqindex
  # y3 <- sagaindex
  y4 <- oemvrindex
  ytwothird <- eqsaga(datasizes)
  ylin <- eqiem(datasizes)
  df <- data.frame(x,y1,y2,y4,ytwothird,ylin)

  print(ggplot(df, aes(x),show.legend = TRUE) +                    
    geom_line(aes(y=y1), colour="red") +
  geom_line(aes(y=y2), colour="blue") +
  geom_line(aes(y=y4), colour="purple") +
  geom_line(aes(y=ytwothird), colour="black") +
  geom_line(aes(y=ylin), colour="brown") +
  scale_y_log10()  + scale_x_log10() +
  xlab("Dataset size") + ylab("Iteration")  +
  ggtitle(precision))

}




epochs
start =2
end = 10

em_ep <- eml
for (i in (1:length(datasizes))){
	em_ep[[i]]$iteration <- datasizes[i]*eml[[i]]$iteration
}


iemseql <- iemseqiternew[[i]][epochs,]
iemseql$iteration <- 1:(K/n)

oemvrl <- oemvriternew[[i]][epochs,]
oemvrl$iteration <- 1:(K/n)


variance <- rbind(eml[[i]][start:end,c(1,4,8)],
                  iemseql[[i]][start:end,c(1,4,8)],
                  sagal[[i]][start:end,c(1,4,8)],
                   oemvrl[[i]][start:end,c(1,4,8)])

graphConvMC2_new(variance, title="IEMs GMM 1e5",legend=TRUE)


length(oemvrl)
curveoemvr <- oemvrl
curvesaga <- sagal
curveem <- eml
curveiem <- ieml
curveiemseq <- iemseql

for (i in (1:length(datasizes))){
	curveoemvr[[i]]$size <- datasizes[i]
  curvesaga[[i]]$size <- datasizes[i]
  curveem[[i]]$size <- datasizes[i]
  curveiem[[i]]$size <- datasizes[i]
  curveiemseq[[i]]$size <- datasizes[i]
}



epochs
start =1
end = 10

plot.oemvr <- NULL
for (i in (1:length(datasizes))){
  plot.oemvr <- rbind(plot.oemvr, curveoemvr[[i]][start:end,c(1,4,9)])

}

plotagainstn(plot.oemvr, title="SAGA-EM for different dataset size",legend=TRUE)


