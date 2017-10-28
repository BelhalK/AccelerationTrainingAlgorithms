library(mlxR)

model2 <- inlineModel("

[LONGITUDINAL]
input = {beta,lambda}  

EQUATION:
h=(beta/lambda)*(t/lambda)^(beta-1)

DEFINITION:
e = {type               = event, 
     rightCensoringTime = 6,  
     hazard             = h}
[INDIVIDUAL]
input={lambda_pop, o_lambda,beta_pop, o_beta}
                      
DEFINITION:
lambda  ={distribution=lognormal, prediction=lambda_pop,  sd=o_lambda}
beta  ={distribution=lognormal, prediction=beta_pop,  sd=o_beta}
     ")


p <- c(lambda_pop=0.3, o_lambda=0.05,
       beta_pop = 1,o_beta = 0.1)
h <- list(name='h', time=seq(0, 6, by=1))
e <- list(name='e', time=0)

N <- 10
res <- simulx(model     = model2, 
              settings  = list(seed=123),
              parameter = p, 
              output    = list(h, e), 
               group     = list(size = N))

print(res$e)



writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent/rtte1.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/timetoevent/rtte1.csv", header=T, sep=","))
