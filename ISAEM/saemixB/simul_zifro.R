library(mlxR)

model <- inlineModel("
              [INDIVIDUAL]
              input = {Tlag_pop, omega_Tlag, ka_pop, omega_ka, V_pop, omega_V, alpha_pop ,omega_alpha, beta_pop, omega_beta}

              DEFINITION:
              Tlag = {distribution=lognormal, typical=Tlag_pop, sd=omega_Tlag}
              ka = {distribution=lognormal, typical=ka_pop, sd=omega_ka}
              V = {distribution=lognormal, typical=V_pop,sd=omega_V}
              alpha = {distribution=lognormal, typical=alpha_pop,sd=omega_alpha}
              beta = {distribution=lognormal, typical=beta_pop,sd=omega_beta}


              [LONGITUDINAL]
              input = {Tlag, ka, V, alpha, beta,a}

              EQUATION:
              Cc = pkmodel(Tlag, ka, V, Cl=alpha*(V^beta))

              OUTPUT:
              output = Cc

              DEFINITION:
              y1 = {distribution=normal, prediction=Cc, errorModel=constant(a)}

                      ")

adm  <- list(amount=100, time=seq(0,50,by=50))


p <- c(Tlag_pop=0.000359, omega_Tlag=3.7,  
      ka_pop=0.7234, omega_ka=0.8,
       V_pop=94.84, omega_V=0.6, 
       alpha_pop=1, omega_alpha=0,  
       beta_pop=0.5, omega_beta=0,  
       a=0.2)
y1 <- list(name='y1', time=seq(1,to=50,by=2))

res <- simulx(model = model,
                 treatment = adm,
                 parameter = p,
                 group = list(size=500, level="individual"),
                 output = y1)
writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/novariability/data/zifro.csv")



# library(mlxR)

# model <- inlineModel("
#               [INDIVIDUAL]
#               input = {Tlag_pop, omega_Tlag, ka_pop, omega_ka, V_pop, omega_V, alpha_pop ,omega_alpha, beta_pop, omega_beta}

#               DEFINITION:
#               Tlag = {distribution=lognormal, typical=Tlag_pop, sd=omega_Tlag}
#               ka = {distribution=lognormal, typical=ka_pop, sd=omega_ka}
#               V = {distribution=lognormal, typical=V_pop,sd=omega_V}
#               alpha = {distribution=lognormal, typical=alpha_pop,sd=omega_alpha}
#               beta = {distribution=lognormal, typical=beta_pop,sd=omega_beta}


#               [LONGITUDINAL]
#               input = {Tlag, ka, V, alpha, beta,a}

#               EQUATION:
#               Cc = pkmodel(Tlag, ka, V, Cl=alpha*(V^beta))

#               OUTPUT:
#               output = Cc

#               DEFINITION:
#               y1 = {distribution=normal, prediction=Cc, errorModel=constant(a)}

#                       ")

# # treatment
# trt <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/zifro/treatment.txt", header = TRUE) 
# # parameters 
# originalId<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/zifro/originalId.txt', header=TRUE) 
# populationParameter<- read.vector('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/zifro/populationParameter_zifro.txt') 
# list.param <- list(populationParameter)
# # output 
# name<-"y1"
# time<-read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/zifro/output1.txt",header=TRUE)
# out1<-list(name=name,time=time) 
# # call the simulator 
# res <- simulx(model=model,treatment=trt,parameter=list.param,output=out1)
# writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/novariability/data/zifro.csv")


# library(mlxR)

# model <- inlineModel("
#             [INDIVIDUAL]
#               input = {Tlag_pop, omega_Tlag, ka_pop, omega_ka, V_pop, omega_V, Cl_pop, omega_Cl}

#               DEFINITION:
#               Tlag = {distribution=lognormal, typical=Tlag_pop, sd=omega_Tlag}
#               ka = {distribution=lognormal, typical=ka_pop, sd=omega_ka}
#               V = {distribution=lognormal, typical=V_pop,sd=omega_V}
#               Cl = {distribution=lognormal, typical=Cl_pop, sd=omega_Cl}


#               [LONGITUDINAL]
#               input =  {Tlag, ka, V, Cl,a}

#               EQUATION:
#               Cc = pkmodel(Tlag, ka, V, Cl)

#               OUTPUT:
#               output = {Cc}

#               DEFINITION:
#               y1 = {distribution=normal, prediction=Cc, errorModel=constant(a)}

#                       ")

# # treatment
# trt <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/zifro/treatment.txt", header = TRUE) 
# # parameters 
# originalId<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/zifro/originalId.txt', header=TRUE) 
# populationParameter<- read.vector('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/zifro/populationParameter2.txt') 
# list.param <- list(populationParameter)
# # output 
# name<-"y1"
# time<-read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/zifro/output1.txt",header=TRUE)
# out1<-list(name=name,time=time) 
# # call the simulator 
# res <- simulx(model=model,treatment=trt,parameter=list.param,output=out1)

# writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/novariability/data/zifro.csv")

