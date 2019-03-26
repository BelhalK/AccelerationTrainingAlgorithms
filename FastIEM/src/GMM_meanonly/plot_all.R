require(ggplot2)
require(gridExtra)
require(reshape2)
library(rlist)

source("utils/algos.R")
source("utils/func.R")
source("utils/plots.R")
theme_set(theme_bw())
options(digits = 22)


# load("RData_VM/precisionagainstn_VM2.RData")
load("RData_VM2/precisionagainstn_VM.RData")



# library("R.utils")
# value <- loadToEnv("RData_VM2/precisionagainstn_VM.RData")[["datasizes"]];

datasizes <- seq(1000, 106000, 5000)
# # save(emiter,file="emiter.RData")
# save(iemiter,file="RData/iemiter.RData")
# save(iemseqiter,file="RData/iemseqiter.RData")
# save(oemvriter,file="RData/oemvriter.RData")
# save(sagaiter,file="RData/sagaiter.RData")
# save(oemiter,file="RData/oemiter.RData")

# iemiternew <- iemiter
# iemseqiternew <- iemseqiter
# oemvriternew <- oemvriter


length(datasizes)
length(oemiter)

em_ep <- eml
for (i in (1:length(datasizes))){
	em_ep[[i]]$iteration <- datasizes[i]*eml[[i]]$iteration
}

eqiem = function(x){x+100}
eqsaga = function(x){x**(2/3)}


for (precision in c(1e-2,1e-3,1e-4,1e-5)){
  emindex  <- iemseqindex <- iemindex  <- oemvrindex <- sagaindex  <- c()

  for (i in (1:length(datasizes))){
  	emindex[i] <- datasizes[[i]]*which(emiter[[i]][,c(4)] < precision)[1]
    iemindex[i] <- which(iemiter[[i]][,c(4)] < precision)[1]
    iemseqindex[i] <- which(iemseqiter[[i]][,c(4)] < precision)[1]
    sagaindex[i] <- which(sagaiter[[i]][,c(4)] < precision)[1]
    oemvrindex[i] <- which(oemvriter[[i]][,c(4)] < precision)[1]
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


