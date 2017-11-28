library(mlxR)

model2 <- inlineModel("
                      [LONGITUDINAL]
                      input = {e, em, ec,a}

                      EQUATION:
                      Cc = e+(em*t)/(ec+t)
                      
                      DEFINITION:
                      y1 ={distribution=normal, prediction=Cc, sd=a}
                      
                      [INDIVIDUAL]
                      input={e_pop,o_e,em_pop,o_em,ec_pop,o_ec}
                      
                      DEFINITION:
                      e   ={distribution=lognormal, prediction=e_pop,   sd=o_e}
                      em   ={distribution=lognormal, prediction=em_pop,   sd=o_em}
                      ec  ={distribution=lognormal, prediction=ec_pop,  sd=o_ec}
                      ")

adm  <- list(amount=c(0, 10, 90), time=c(0, 10, 90))


p <- c(e_pop=20, o_e=0.5,
       em_pop=100, o_em=0.5, 
       ec_pop=10, o_ec=0.5,  
       a=0.1)
y1 <- list(name='y1', time=c(0, 10, 90))


res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=50, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/pd1.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/pd1.csv", header=T, sep=","))


#modification for mlxsaem dataread function
obj <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/pd1.csv", header=T, sep=",")
a <- which(obj[,4]==0)
data <- obj[-a,]
b <- which(data[,4]==10)
data <- data[-b,]
c <- which(data[,4]==90)
data <- data[-c,]

write.table(data, "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/pd1.csv", sep=",", row.names=FALSE,quote = FALSE, col.names=TRUE)
