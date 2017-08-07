#library(rstan)
setwd("/Users/karimimohammedbelhal/Desktop/variationalBayes/mcmc_R_isolate/Dir2")
  source('compute_LL.R') 
  source('func_aux.R') 
  source('func_cov.R') 



index = 1
graphConvMC_twokernels(post_rwm[[index]],post_rwm[[index]], title="rwm vs foce")
graphConvMC_twokernels(post_rwm[[index]],post_foce[[index]], title="rwm vs foce")


final_rwm <- post_rwm[[1]]
for (i in 2:length(post_rwm)) {
  final_rwm <- rbind(final_rwm, post_rwm[[i]])
}


final_foce <- post_foce[[1]]
for (i in 2:length(post_foce)) {
  final_foce <- rbind(final_foce, post_foce[[i]])
}



graphConvMC_twokernels(final_rwm,final_rwm, title="EM")
graphConvMC_twokernels(final_rwm,final_foce, title="EM")


#Autocorrelation
rwm.obj <- as.mcmc(post_rwm[[1]])
corr_rwm <- autocorr(rwm.obj[,2])
autocorr.plot(rwm.obj[,2])

foce.obj <- as.mcmc(post_foce[[1]])
corr_foce <- autocorr(foce.obj[,2])
autocorr.plot(foce.obj[,2])


#MSJD
mssd(post_rwm[[index]][,2])
mssd(post_foce[[index]][,2])



