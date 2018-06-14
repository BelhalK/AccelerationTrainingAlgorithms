# File generated automatically on 2017-10-31 12:50:43
 
library(mlxR)  
 
setwd(dirname(parent.frame(2)$ofile)) 

# model 
model<-"rtteWeibull_project_model.txt"

# parameters 
originalId<- read.table('originalId.txt', header=TRUE) 
populationParameter<- read.vector('populationParameter.txt') 
list.param <- list(populationParameter)
# output 
name<-"Event"
time<-read.table("output.txt",header=TRUE)
out<-list(name=name,time=time) 

# call the simulator 
res <- simulx(model=model,parameter=list.param,output=out)
