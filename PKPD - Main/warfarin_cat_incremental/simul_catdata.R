library(mlxR)

model2 <- inlineModel("

              [LONGITUDINAL]
              input = {th1, th2, th3}

              EQUATION:
              lgp0 = th1
              lgp1 = lgp0 + th2
              lgp2 = lgp1 + th3

              DEFINITION:
              level = { type = categorical,  categories = {0, 1, 2, 3},
              logit(P(level<=0)) = th1
              logit(P(level<=1)) = th1 + th2
              logit(P(level<=2)) = th1 + th2 + th3
              }

              [INDIVIDUAL]
              input={th1_pop, o_th1,th2_pop, o_th2,th3_pop, o_th3}
                      

              DEFINITION:
              th1  ={distribution=normal, prediction=th1_pop,  sd=o_th1}
              th2  ={distribution=lognormal, prediction=th2_pop,  sd=o_th2}
              th3  ={distribution=lognormal, prediction=th3_pop,  sd=o_th3}
                      
                      ")

p <- c(th1_pop=1, o_th1=0.5,
       th2_pop=2, o_th2=0.8, 
       th3_pop=2.5, o_th3=0.6)



y1 <- list(name='level', time=seq(1,to=100,by=10))



res2a2 <- simulx(model = model2,
                 parameter = p,
                 group = list(size=100, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat1.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/warfarin_cat/data/cat1.csv", header=T, sep=","))


