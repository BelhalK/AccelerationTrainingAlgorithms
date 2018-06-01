############################### Simulation - MCMC kernels (E-step) #############################

indiv.variational.inference<-function(model,data,control=list()) {
	# E-step - simulate unknown parameters
	# Input: kiter, Uargs, structural.model, mean.phi (unchanged)
	# Output: varList, DYF, phiM (changed)
	kiter <- 1
	saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
	saemix.options<-saemixObject["options"]
	saemix.model<-saemixObject["model"]
	# Initialising random generator
	OLDRAND<-TRUE
	set.seed(saemix.options$seed)
	#intitialisation
	# xinit<-initialiseMainAlgo(saemix.data,saemix.model,saemix.options)
	xinit<-initialiseMainAlgo(saemix.data,saemix.model,saemix.options)
	saemix.model<-xinit$saemix.model
	Dargs<-xinit$Dargs
	Uargs<-xinit$Uargs
	varList<-xinit$varList
	phiM<-xinit$phiM
	mean.phi<-xinit$mean.phi
	DYF<-xinit$DYF
	opt<-xinit$opt
	structural.model<-saemix.model["model"]

	# Function to perform MCMC simulation
	nb.etas<-length(varList$ind.eta)
	domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
	omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
	omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
	somega <- solve(omega.eta)
	# "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)

	VK<-rep(c(1:nb.etas),2)
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
	saemix.options<-saemixObject["options"]

	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	phiMc<-phiM

	#Initialization
	L <- 10 #nb iterations MONTE CARLO
	rho <- 10e-12 #gradient ascent stepsize
	K <- control$nbiter.gd

	#VI to find the right mean mu (gradient descent along the elbo)
	i <- 10
	mu <- list(etaM[i,],etaM[i,])
	
	mu[[1]][1] = 20
	mu[[1]][2] = 0.1
	#if Gamma fixed to Laplace Gamma
	trueGamma <- control$Gamma.laplace
	chol.Gamma <- chol(trueGamma[[i]])
	inv.Gamma <- solve(trueGamma[[i]])

	for (k in 1:K) {
		if (k%%10==0) print(k)
			sample <- list(etaM[i,],etaM[i,])  #list of samples for monte carlo integration
			sample1 <- list(etaM[i,],etaM[i,])  #list of samples for gradient computation
			estim <- list(etaM[i,],etaM[i,])
			gradlogq <- etaM[i,]
			for (l in 1:L) {
				sample[[l]] <- mu[[k]] +matrix(rnorm(nb.etas), ncol=nb.etas)%*%chol.Gamma
				phiMc[i,varList$ind.eta]<-mean.phiM[i,varList$ind.eta]+sample[[l]]
				psiMc<-transphi(phiMc,Dargs$transform.par)
				fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred<-log(cutoff(fpred))
				gpred<-error(fpred,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
				#Log complete computation
				sumDYF <- colSums(DYF)
				logp <- sumDYF[i] + 0.5*rowSums(sample[[l]]*(sample[[l]]%*%somega))
				#Log proposal computation
				logq <- 0.5*rowSums(sample[[l]]*(sample[[l]]%*%inv.Gamma))
				#gradlogq computation
				for (j in 1:nb.etas) {
					mu2 <- mu[[k]]
					mu2[j] <- mu[[k]][j] + mu[[k]][j]/100
					sample2 <- mu2 +matrix(rnorm(nb.etas), ncol=nb.etas)%*%chol.Gamma
					gradlogq[j] <- (0.5*rowSums(sample2*(sample2%*%inv.Gamma)) - 0.5*rowSums(sample[[l]]*(sample[[l]]%*%inv.Gamma))) / (mu[[k]][j]/100)
				}
				estim[[l]] <- (logq - logp)*gradlogq
			}
			grad_mu_elbo <- 1/L*Reduce("+", estim) 
			#Gradient ascent along that gradient
			mu[[k+1]] <- mu[[k]] - rho*grad_mu_elbo
	}

	mu.vi <- mu[[K]]
	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM[,varList$ind.eta]
	return(list(mu=mu))
}
