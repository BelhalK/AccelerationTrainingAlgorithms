library(mlxR)

model2 <- inlineModel("
                      [LONGITUDINAL]
                      input = {ymax, xmax, slope, a}
                      
                      EQUATION:
                      Cc = ymax+slope*(t-xmax)
                      
                      DEFINITION:
                      y1 ={distribution=normal, prediction=Cc, sd=a}
                      
                      [INDIVIDUAL]
                      input={ymax_pop, o_ymax,xmax_pop, o_xmax,slope_pop, o_slope}
                      
                      DEFINITION:
                      ymax  ={distribution=normal, prediction=ymax_pop,  sd=o_ymax}
                      xmax  ={distribution=normal, prediction=xmax_pop,  sd=o_xmax}
                      slope  ={distribution=normal, prediction=slope_pop,  sd=o_slope}
                      
                      ")

adm  <- list(amount=1, time=seq(0,50,by=50))
p <- c(ymax_pop=10, o_ymax=0.5,
       xmax_pop=1, o_xmax=1, 
       slope_pop=2, o_slope=0.1,
       a=0.1)
y1 <- list(name='y1', time=seq(1,to=50,by=10))


res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=20, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/yield/yield_synth.csv")
table <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/yield/yield_synth.csv", header=T, sep=",")
head(table)
table[1:45,]


#modification for mlxsaem dataread function
obj <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/yield/yield_synth.csv", header=T, sep=",")
obj <- obj[obj$amount !=1,]
obj[,4] <- 1
write.table(obj, "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/yield/yield_synth.csv", sep=",", row.names=FALSE,quote = FALSE, col.names=TRUE)