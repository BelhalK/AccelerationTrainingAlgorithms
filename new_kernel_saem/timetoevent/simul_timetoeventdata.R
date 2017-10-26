library(mlxR)

model2 <- inlineModel("

[LONGITUDINAL]
input = {beta,lambda}  

EQUATION:
h=lambda/beta

DEFINITION:
e = {type               = event, 
     maxEventNumber     = 5, 
     rightCensoringTime = 60,  
     hazard             = h}
     ")


p <- c(beta = 1, lambda=20)
h <- list(name='h', time=seq(0, 6, by=1))
e <- list(name='e', time=0)

N <- 50
res <- simulx(model     = model2, 
              settings  = list(seed=123),
              parameter = p, 
              output    = list(h, e), 
               group     = list(size = N))

print(res$e)



writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent/timeto.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent/timeto.csv", header=T, sep=","))
