library(FactoMineR)
ozone <- read.table("data/ozone.txt", header = T, sep = " ", row.names = 1)
res.pca <- PCA(ozone, quanti.sup = 1, quali.sup = c(12, 13))