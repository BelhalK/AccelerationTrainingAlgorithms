library(mlxR)


model2 <- inlineModel("
                      [LONGITUDINAL]
                      input = {ka, V, k, a}
                      EQUATION:
                      Cc = ka/(V*(ka-k))*(exp(-k*t)-exp(-ka*t))
                      
                      DEFINITION:
                      y1 ={distribution=normal, prediction=Cc, sd=a}
                      
                      [INDIVIDUAL]
                      input={ka_pop,o_ka, V_pop,o_V, k_pop,o_k}
                      
                      DEFINITION:
                      ka  ={distribution=lognormal, prediction=ka_pop,  sd=o_ka}
                      V  ={distribution=lognormal, prediction=V_pop,  sd=o_V}
                      k  ={distribution=lognormal, prediction=k_pop,  sd=o_k}                      
                      ")
adm  <- list(amount=1000, time=seq(0,50,by=50))
p <- c(ka_pop=1, o_ka=0.5,
       V_pop=20, o_V=0.2, 
       k_pop=2, o_k=0.1,  
       a=0.1)
y1 <- list(name='y1', time=seq(1,to=50,by=5))


res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=10, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/theo/theo_synth.csv")
table <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/theo/theo_synth.csv", header=T, sep=",")
head(table)
table[1:45,]


#modification for mlxsaem dataread function
obj <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/theo/theo_synth.csv", header=T, sep=",")
obj <- obj[obj$amount !=1,]
obj[,4] <- 1000
write.table(obj, "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/theo/theo_synth.csv", sep=",", row.names=FALSE,quote = FALSE, col.names=TRUE)


