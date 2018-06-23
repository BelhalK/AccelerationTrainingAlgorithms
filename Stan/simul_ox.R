library("mlxR")
library(saemix)

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

model <- inlineModel("


[INDIVIDUAL]
input = {base_pop, slope_pop,  omega_base, omega_slope}
DEFINITION:
base = {distribution=normal, reference=base_pop, sd=omega_base}
slope  = {distribution=normal, reference=slope_pop,  sd=omega_slope }


[LONGITUDINAL]
input = {base, slope,a}
EQUATION:
C = base + slope*t
DEFINITION:
y = {distribution=normal, prediction=C, sd=a}
")

N=100

param   <- c(
  base_pop  = 140,    omega_base  = 0.5,
  slope_pop   = 1,   omega_slope   = 0.4, a =1)
  
res <- simulx(model     = model,
              parameter = param,
              group     = list(size=N, level='individual'),
              output    = list(name='y', time=seq(1,3,by=1)))

 head(res$y)


