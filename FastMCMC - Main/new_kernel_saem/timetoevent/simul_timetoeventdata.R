library(mlxR)

model2 <- inlineModel("

[LONGITUDINAL]
input = {beta,lambda}  

EQUATION:
h=lambda/beta

DEFINITION:
e = {type               = event, 
     rightCensoringTime = 60,  
     hazard             = h}
[INDIVIDUAL]
input={lambda_pop, o_lambda}
                      
DEFINITION:
lambda  ={distribution=lognormal, prediction=lambda_pop,  sd=o_lambda}
     ")


p <- c(lambda_pop=0.3, o_lambda=0.5,
       beta = 1)
h <- list(name='h', time=seq(0, 60, by=1))
e <- list(name='e', time=0)

N <- 10
res <- simulx(model     = model2, 
              settings  = list(seed=123),
              parameter = p, 
              output    = list(h, e), 
               group     = list(size = N))

print(res$e)



writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent/timeto.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent/timeto.csv", header=T, sep=","))
