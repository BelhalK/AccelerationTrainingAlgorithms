library(mlxR)


model2 <- inlineModel("
                      [LONGITUDINAL]
                      input = {ka, V, k, a}
                      EQUATION:
                      Cc = ka/(V*(ka-k))*(exp(-k*t)-exp(-ka*t))
                      
                      DEFINITION:
                      y1 ={distribution=lognormal, prediction=Cc, sd=a}
                      
                      [INDIVIDUAL]
                      input={ka_pop,o_ka, V_pop,o_V, k_pop,o_k}
                      
                      DEFINITION:
                      ka  ={distribution=lognormal, prediction=ka_pop,  sd=o_ka}
                      V  ={distribution=lognormal, prediction=V_pop,  sd=o_V}
                      k  ={distribution=lognormal, prediction=k_pop,  sd=o_k}                      
                      ")

adm  <- list(amount=300, time=seq(0,100,by=1))
p <- c(ka_pop=1, o_ka=0.5,
       V_pop=20, o_V=1, 
       k_pop=2, o_k=0.1,  
       a=0.1)
y1 <- list(name='y1', time=seq(1,to=100,by=10))


res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=10, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/theo/theo_synth.csv")
head(read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/theo/theo_synth.csv", header=T, sep=","))

#modification for mlxsaem dataread function
obj <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/mcmc_newkernel_saem/theo/theo_synth.csv", header=T, sep=";")
obj <- obj[obj$amount !=1,]
write.table(obj, "/Users/karimimohammedbelhal/Documents/GitHub/saem/mcmc_newkernel_saem/theo/theonew.csv", sep=",", row.names=FALSE,quote = FALSE, col.names=TRUE)


growth.linear<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (2 columns, base and slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  base<-psi[id,1]
  slope<-psi[id,2]
  f<-base+slope*x
  return(f)
}
saemix.model<-saemixModel(model=growth.linear,description="Linear model",
  psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(1,0),covariance.model=matrix(c(1,1,1,1),ncol=2,byrow=TRUE), 
  error.model="constant")