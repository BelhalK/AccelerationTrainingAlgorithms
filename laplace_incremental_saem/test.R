############################### Simulation - MCMC kernels (E-step) #############################
estep_laplace_incremental<-function(kiter, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM,saemixObject) {
	# E-step - simulate unknown parameters
	# Input: kiter, Uargs, structural.model, mean.phi (unchanged)
	# Output: varList, DYF, phiM (changed)
	
	# Function to perform MCMC simulation

	nb.etas<-length(varList$ind.eta)
	domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
	omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
	omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
	chol.omega<-try(chol(omega.eta))
	somega<-solve(omega.eta)
	saemix.options<-saemixObject["options"]
	# "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)
	
	VK<-rep(c(1:nb.etas),2)
	Uargs$nchains = 1
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
	psiM<-transphi(phiM,Dargs$transform.par)
	fpred<-structural.model(psiM, Dargs$IdM, Dargs$X
	phiMc<-phiM
	# map_range <- c(1:3)
	map_range <- saemix.options$map.range
	# map_range <- c(15,25,	phiM[ind_rand,varList$ind.eta] <- phiMold[ind_rand,varList$ind.eta]

	# phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM, post = post))
}
