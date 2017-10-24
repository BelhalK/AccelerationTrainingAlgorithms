library(mlxR)

model2 <- inlineModel("

[LONGITUDINAL]
input = {lambda,beta}  

EQUATION:
h=lambda*beta

DEFINITION:
e = {type               = event, 
     eventType          = intervalCensored, 
     intervalLength     = 5, 
     rightCensoringTime = 60,  
     hazard             = h}



                      ")


p <- c(lambda=20, beta=1)
h <- list(name='h', time=seq(0, 60, by=10))
e <- list(name='e', time=0)

N <- 10
res2a2 <- simulx(model     = model2, 
              settings  = list(seed=1234),
              parameter = p, 
              output    = list(h, e), 
              group     = list(size = N))

writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent/timeto.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent/timeto.csv", header=T, sep=","))


