library(MASS)
library(ggplot2)
library(reshape2)

load("RData/trauma.RData")

#PLOTS
dim = 1
imcem = list.imcem$mu[dim,]
saem = list.saem$mu[dim,]
mcem = list.mcem$mu[dim,]

# imcem = list.imcem$seqbeta[dim,]
# saem = list.saem$seqbeta[dim,]
# mcem = list.mcem$seqbeta[dim,]

x = 1:length(imcem)
df <- data.frame(x,saem,mcem, imcem)
ggplot(data=df)+
  geom_line(mapping=aes(y=saem,x= x,color="saem"),size=0.5 ) +
  geom_line(mapping=aes(y=mcem,x= x,color="mcem"),size=0.5) +
  geom_line(mapping=aes(y=imcem,x= x,color="imcem"),size=0.5) +
  scale_color_manual(values = c(
    'saem' = 'darkblue','mcem' = 'red','imcem' = 'black')) +
  labs(color = 'Algo')+ ylab("beta")


#PER EPOCHS
epochs = seq(1,nb.iter,N/batchsize)
x = 2:(nb.epochs+1)
saem.ep <- saem[x]
mcem.ep <- mcem[x]
imcem.ep <- imcem[(epochs+1)]
df <- data.frame(x,saem.ep,mcem.ep,imcem.ep)

ggplot(data=df)+
  geom_line(mapping=aes(y=saem.ep,x= x,color="saem"),size=0.5 ) +
  geom_line(mapping=aes(y=mcem.ep,x= x,color="mcem"),size=0.5) +
  geom_line(mapping=aes(y=imcem.ep,x= x,color="imcem"),size=0.5) +
  scale_color_manual(values = c(
    'saem' = 'darkblue','mcem' = 'red','imcem' = 'black')) +
  labs(color = 'Algo')