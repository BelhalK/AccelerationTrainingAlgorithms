library("mlxR")

 
Tlag_true=0.78
ka_true <- 1
V_true <- 8
Cl_true <- 0.1

o_Tlag <- 0.57 #o^2=0.32
o_ka <- 0.5 #o^2=0.25
o_V <- 0.2  #o^2=0.04
o_Cl <- 0.3  #o^2=0.09
a_true = 0.266
beta_Cl_lw70_true = 0.60411
beta_V_lw70_true = 0.8818
seed0 = 39546

  myModel <- inlineModel("

[COVARIATE]
input = wt

EQUATION:
lw70 = log(wt/70)

[INDIVIDUAL]
input = {Tlag_pop, omega_Tlag, ka_pop, omega_ka, V_pop, beta_V_lw70, lw70, omega_V, Cl_pop, beta_Cl_lw70, omega_Cl}

DEFINITION:
Tlag = {distribution=lognormal, typical=Tlag_pop, sd=omega_Tlag}
ka = {distribution=lognormal, typical=ka_pop, sd=omega_ka}
V = {distribution=lognormal, typical=V_pop, covariate=lw70, coefficient=beta_V_lw70, sd=omega_V}
Cl = {distribution=lognormal, typical=Cl_pop, covariate=lw70, coefficient=beta_Cl_lw70, sd=omega_Cl}

[LONGITUDINAL]
input =  {Tlag, ka, V, Cl,a}

EQUATION:
Cc = pkmodel(Tlag, ka, V, Cl)

OUTPUT:
output = {Cc}

DEFINITION:
y1 = {distribution=normal, prediction=Cc, sd=a}
")


populationParameter   <- c(Tlag_pop= Tlag_true, omega_Tlag= o_Tlag,
  ka_pop  = ka_true,    omega_ka  = o_ka,
  V_pop   = V_true,   omega_V   = o_V,
  Cl_pop  = Cl_true,    omega_Cl  = o_Cl, a =a_true, beta_V_lw70 = beta_V_lw70_true, beta_Cl_lw70 = beta_Cl_lw70_true)

trt <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/design2/treatment.txt", header = TRUE) 
originalId<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/design2/originalId.txt', header=TRUE) 
individualCovariate<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/design2/individualCovariate.txt', header = TRUE) 
time<-read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/design2/output1.txt",header=TRUE)

# trt <- read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/treatment.txt", header = TRUE) 
# originalId<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/originalId.txt', header=TRUE) 
# individualCovariate<- read.table('/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/individualCovariate.txt', header = TRUE) 
# time<-read.table("/Users/karimimohammedbelhal/Desktop/CSDA_code_ref/warfarin/output1.txt",header=TRUE)


list.param <- list(populationParameter,individualCovariate)
name<-"y1"
out1<-list(name=name,time=time) 

# call the simulator 
res <- simulx(model=myModel,treatment=trt,parameter=list.param,output=out1)
# writeDatamlx(res, result.file = paste("/Users/karimimohammedbelhal/Desktop/data_pk/pk_mcstudy_", m, ".csv", sep=""))
# warfarin.saemix<-read.table(paste("/Users/karimimohammedbelhal/Desktop/data_pk/pk_mcstudy_", m, ".csv", sep=""), header=T, sep=",")
# typeof(warfarin.saemix[2,3])
# typeof(res$y1[2,3])
individualCovariate$wt <- log(individualCovariate$wt/70)
warfarin.saemix <- res$y1
treat <- res$treatment[,c(1,3)]
covandtreat <- merge(individualCovariate ,treat,by="id")
warfarin.saemix <- merge(covandtreat ,warfarin.saemix,by="id")


# writeDatamlx(res, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/pkcov.csv")



# individualCovariate$wt <- log(individualCovariate$wt/70)
# warfarin.saemix <- res$y1
# treat <- res$treatment[,c(1,3)]
# covandtreat <- merge(individualCovariate ,treat,by="id")
# warfarin.saemix <- merge(covandtreat ,warfarin.saemix,by="id")



