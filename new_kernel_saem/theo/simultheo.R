library(mlxR)

model2 <- inlineModel("
                      [LONGITUDINAL]
                      input = {A1, A2, A3, alpha1, alpha2, alpha3, a}
                      
                      EQUATION:
                      Cc = A1*exp(-alpha1*t)+A2*exp(-alpha2*t)+A3*exp(-alpha3*t)
                      
                      DEFINITION:
                      y1 ={distribution=lognormal, prediction=Cc, sd=a}
                      
                      [INDIVIDUAL]
                      input={A1_pop, o_A1,A2_pop, o_A2,A3_pop, o_A3,alpha1_pop, o_alpha1,alpha2_pop, o_alpha2,alpha3_pop, o_alpha3}
                      
                      DEFINITION:
                      A1  ={distribution=lognormal, prediction=A1_pop,  sd=o_A1}
                      A2  ={distribution=lognormal, prediction=A2_pop,  sd=o_A2}
                      A3  ={distribution=lognormal, prediction=A3_pop,  sd=o_A3}
                      alpha1  ={distribution=lognormal, prediction=alpha1_pop,  sd=o_alpha1}
                      alpha2  ={distribution=lognormal, prediction=alpha2_pop,  sd=o_alpha2}
                      alpha3  ={distribution=lognormal, prediction=alpha3_pop,  sd=o_alpha3}
                      
                      ")

adm  <- list(amount=1, time=seq(0,50,by=50))
p <- c(A1_pop=60, o_A1=0.5,
       A2_pop=1, o_A2=1, 
       A3_pop=100, o_A3=0.1,  
       alpha1_pop=0.6,  o_alpha1=0,
       alpha2_pop=10,  o_alpha2=0.3,
       alpha3_pop=0.1,  o_alpha3=0,
       a=0.1)
y1 <- list(name='y1', time=seq(1,to=50,by=2))


res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=500, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/mcmc_newkernel_saem/theo/theo.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/mcmc_newkernel_saem/theo/theo.csv", header=T, sep=","))

#modification for mlxsaem dataread function
obj <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/mcmc_newkernel_saem/theo/theo.csv", header=T, sep=";")
obj <- obj[obj$amount !=1,]
write.table(obj, "/Users/karimimohammedbelhal/Documents/GitHub/saem/mcmc_newkernel_saem/theo/theonew.csv", sep=",", row.names=FALSE,quote = FALSE, col.names=TRUE)
