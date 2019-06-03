library(MASS)
library(ggplot2)
library(reshape2)

#load("RData/trauma_all.RData")
load("RData/trauma_final.RData")
#PLOTS
dim = 3
imcem10 = list.imcem10$mu[dim,]
imcem50 = list.imcem50$mu[dim,]
imcem = list.imcem$mu[dim,]
saem = list.saem$mu[dim,]
mcem = list.mcem$mu[dim,]

plot(saem)

plot(imcem)
plot(mcem)
plot(imcem10)
plot(saem)
# imcem = list.imcem$seqbeta[dim,]
# saem = list.saem$seqbeta[dim,]
# mcem = list.mcem$seqbeta[dim,]

x = 1:length(saem)
plot(x,mcem)
df <- data.frame(x,saem,mcem)
ggplot(data=df)+
  geom_line(mapping=aes(y=saem,x= x,color="saem"),size=0.5 ) +
  geom_line(mapping=aes(y=mcem,x= x,color="mcem"),size=0.5) +
  scale_color_manual(values = c(
    'saem' = 'darkblue','mcem' = 'red')) +
  labs(color = 'Algo')+ ylab("beta")



x = 1:length(saem)
plot(x,mcem)
df <- data.frame(x,saem,mcem, imcem[x])
ggplot(data=df)+
  geom_line(mapping=aes(y=saem,x= x,color="saem"),size=0.5 ) +
  geom_line(mapping=aes(y=mcem,x= x,color="mcem"),size=0.5) +
  geom_line(mapping=aes(y=imcem,x= x,color="imcem"),size=0.5) +
  scale_color_manual(values = c(
    'saem' = 'darkblue','mcem' = 'red','imcem' = 'black')) +
  labs(color = 'Algo')+ ylab("beta")



#PLOTS MU and BETA


dim = 1
imcem10 = list.imcem10$seqbeta[dim,]
imcem50 = list.imcem50$seqbeta[dim,]
imcem = list.imcem$seqbeta[dim,]
saem = list.saem$seqbeta[dim,]
mcem = list.mcem$seqbeta[dim,]


dim=3
imcem10 = list.imcem10$mu[dim,]
imcem50 = list.imcem50$mu[dim,]
imcem = list.imcem$mu[dim,]
saem = list.saem$mu[dim,]
mcem = list.mcem$mu[dim,]
plot(saem)


#PER EPOCHS
epochs = seq(1,nb.iter,N/batchsize)
nb.epochs = length(epochs)
x = 2:(nb.epochs+1)
saem.ep <- saem[x]
mcem.ep <- mcem[x]
imcem.ep <- imcem[(epochs+1)]


# df <- data.frame(x,saem.ep,mcem.ep,imcem.ep)

# ggplot(data=df)+
#   geom_line(mapping=aes(y=saem.ep,x= x,color="saem"),size=0.5 ) +
#   geom_line(mapping=aes(y=mcem.ep,x= x,color="mcem"),size=0.5) +
#   geom_line(mapping=aes(y=imcem.ep,x= x,color="imcem"),size=0.5) +
#   scale_color_manual(values = c(
#     'saem' = 'darkblue','mcem' = 'red','imcem' = 'black')) +
#   labs(color = 'Algo')


epochs10 = seq(1,10*nb.epochs,10)
epochs50 = seq(1,2*nb.epochs,2)
imcem10.ep <- imcem10[(epochs10+1)]
imcem50.ep <- imcem50[(epochs50+1)]

df <- data.frame(x,saem.ep,mcem.ep,imcem.ep, imcem10.ep, imcem50.ep)

ggplot(data=df)+
  geom_line(mapping=aes(y=saem.ep,x= x,color="saem"),size=0.5 ) +
  geom_line(mapping=aes(y=mcem.ep,x= x,color="mcem"),size=0.5) +
  geom_line(mapping=aes(y=imcem.ep,x= x,color="imcem"),size=0.5) +
  geom_line(mapping=aes(y=imcem10.ep,x= x,color="imcem10"),size=0.5) +
  geom_line(mapping=aes(y=imcem50.ep,x= x,color="imcem50"),size=0.5) +
  scale_color_manual(values = c(
    'saem' = 'darkblue','mcem' = 'red','imcem' = 'black','imcem10' = 'pink','imcem50' = 'green')) +
  labs(color = 'Algo')


write.csv(df, file = "notebooks/mu.csv")

