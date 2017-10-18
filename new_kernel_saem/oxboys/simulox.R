library(mlxR)


model2 <- inlineModel("
                      [LONGITUDINAL]
                      input = {base, slope,a}
                      EQUATION:
                      Cc = base+slope*t

                      DEFINITION:
                      y1 ={distribution=lognormal, prediction=Cc, sd=a}
                      
                      [INDIVIDUAL]
                      input={base_pop,o_base, slope_pop,o_slope}
                      
                      DEFINITION:
                      base  ={distribution=lognormal, prediction=base_pop,  sd=o_base}
                      slope  ={distribution=normal, prediction=slope_pop,  sd=o_slope}
                      ")

adm  <- list(amount=1, time=seq(0,50,by=50))
p <- c(base_pop=100, o_base=2,
       slope_pop=1, o_slope=1,  
       a=1)
y1 <- list(name='y1', time=seq(1,to=50,by=5))


res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=10, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/oxboys/ox_synth.csv")
table <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/oxboys/ox_synth.csv", header=T, sep=",")
head(table)
table[1:45,]


#modification for mlxsaem dataread function
obj <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/oxboys/ox_synth.csv", header=T, sep=",")
obj <- obj[obj$amount !=1,]
obj[,4] <- 1
write.table(obj, "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/oxboys/ox_synth.csv", sep=",", row.names=FALSE,quote = FALSE, col.names=TRUE)


