############################### Simulation - MCMC kernels (E-step) #############################

variational.inference<-function(model,data,control=list()) {
	# E-step - simulate unknown parameters
	# Input: kiter, Uargs, structural.model, mean.phi (unchanged)
	# Output: varList, DYF, phiM (changed)
  kiter <- 1

  saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
  saemix.options<-saemixObject["options"]
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
  saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
  saemix.data@ntot.obs<-dim(saemix.data@data)[1]
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
 	betas<-betas.ini<-xinit$betas
  	fixed.psi<-xinit$fixedpsi.ini
  	var.eta<-varList$diag.omega
  	structural.model<-saemix.model["model"]

	# Function to perform MCMC simulation
	nb.etas<-length(varList$ind.eta)
	domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
	omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
	omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
	chol.omega<-try(chol(omega.eta))
	somega<-solve(omega.eta)
	
	# "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)
	
	VK<-rep(c(1:nb.etas),2)
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
	saemix.options<-saemixObject["options"]
	map_range <- saemix.options$map.range

	if(Dargs$type=="structural"){
		U.y<-compute.LLy_c(phiM,varList$pres,Uargs,Dargs,DYF)
	} else{
		U.y <- compute.LLy_d(phiM,Uargs,Dargs,DYF)
	}
	
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	phiMc<-phiM

U.eta<-0.5*rowSums(etaM*(etaM%*%somega))

#VI with linear model
propc <- U.eta
prop <- U.eta
nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
nrs2<-1

#Initialization

L <- 2 #nb iterations MONTE CARLO
rho <- 0.00000000001 #gradient ascent stepsize
#VI to find the right mean mu (gradient descent along the elbo)

mu <- list(etaM,etaM)
outputGamma <- list(omega.eta,omega.eta)
for (i in 1:Dargs$NM) {
	outputGamma[[i]] <- list(omega.eta,omega.eta)
}
K <- control$nbiter.gd

for (k in 1:K) {
	print(k)
	sample <- list(etaM,etaM)  #list of samples for monte carlo integration
	sample1 <- list(etaM,etaM)  #list of samples for gradient computation
	estim <- list(etaM,etaM)
	gradlogq <- etaM
	for (i in 1:Dargs$NM) {
		sGamma <- outputGamma[[i]][[k]]
		for (l in 1:L) {
			sample[[l]] <- mu[[k]] +matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)%*%chol(outputGamma[[i]][[k]])
			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+sample[[l]]
			psiMc<-transphi(phiMc,Dargs$transform.par)
			fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred<-log(cutoff(fpred))
			gpred<-error(fpred,varList$pres)
			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
			#Log complete computation
			logp <- colSums(DYF) + 0.5*rowSums(sample[[l]]*(sample[[l]]%*%somega))
			#Log proposal computation
			logq <- 0.5*rowSums(sample[[l]]*(sample[[l]]%*%sGamma))
			#gradlogq computation
			for (j in 1:nb.etas) {
				sample1[[l]] <- sample[[l]]
				sample1[[l]][,j] <- sample[[l]][,j] + 0.01
				gradlogq[,j] <- (0.5*rowSums(sample1[[l]]*(sample1[[l]]%*%sGamma)) - 0.5*rowSums(sample[[l]]*(sample[[l]]%*%sGamma))) / 0.01
			}
			estim[[l]] <- sample[[l]]
			estim[[l]][i,] <- (logp[i] - logq[i])*gradlogq[i,]
		}
	}
	grad_mu_elbo <- 1/L*Reduce("+", estim) 
	#Gradient ascent along that gradient
	mu[[k+1]] <- mu[[k]] + rho*grad_mu_elbo
	
	#Update the proposal covariance
	sample <- list(etaM,etaM)  #list of samples for monte carlo integration
	sample1 <- list(etaM,etaM)  #list of samples for gradient computation
	gradlogq <- list(omega.eta,omega.eta)
	for (i in 1:Dargs$NM) {
		gradlogq[[i]] <- omega.eta
	}
	estimcov <- gradlogq

	for (i in 1:Dargs$NM) {
		estimcov[[i]] <- list(omega.eta,omega.eta)
	}

	for (l in 1:L) {
		sample[[l]] <- mu[[k]] +matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)%*%chol(outputGamma[[i]][[k]])
		for (r in 1:nb.etas) {
			Gamma1 <- outputGamma[[i]][[k]]
			for (j in 1:nb.etas) {
				Gamma1[r,] <- outputGamma[[i]][[k]][r,]
				Gamma1[r,j] <- Gamma1[r,j] + 0.01
				sGamma1 <- solve(Gamma1)
				sample1[[l]] <- mu[[k]] +matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)%*%chol(Gamma1)
				
				temp <- (0.5*rowSums(sample1[[l]]*(sample1[[l]]%*%sGamma1)) - 0.5*rowSums(sample[[l]]*(sample[[l]]%*%sGamma1))) / 0.01
				for (i in 1:Dargs$NM) {
					gradlogq[[i]][r,j] <- temp[i]
				}
			}
		}

		phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+sample[[l]]
		psiMc<-transphi(phiMc,Dargs$transform.par)
		fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
		if(Dargs$error.model=="exponential")
			fpred<-log(cutoff(fpred))
		gpred<-error(fpred,varList$pres)
		DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
		#Log complete computation
		logp <- colSums(DYF) + 0.5*rowSums(sample[[l]]*(sample[[l]]%*%somega))
		#Log proposal computation
		logq <- 0.5*rowSums(sample[[l]]*(sample[[l]]%*%sGamma))
		
		
		for (i in 1:Dargs$NM) {
			estimcov[[i]][[l]] <- (logp[i] - logq[i])*gradlogq[[i]]
		}
	}
	grad_cov_elbo <- gradlogq
	for (i in 1:Dargs$NM) {
		grad_cov_elbo[[i]] <- 1/L*Reduce("+", estimcov[[i]])
		outputGamma[[i]][[k+1]] <- outputGamma[[i]][[k]] + rho*grad_cov_elbo[[i]]
	}

}

mu.vi <- mu[[K]]
Gamma.vi <- chol.Gamma.vi <- inv.Gamma.vi <- list(omega.eta,omega.eta)
for (i in 1:(Dargs$NM)){
	Gamma.vi[[i]] <- outputGamma[[i]][[K]]
	# chol.Gamma.vi[[i]] <- chol(Gamma_.vi[[i]])
	# inv.Gamma.vi[[i]] <- solve(Gamma_.vi[[i]])
}

	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM[,varList$ind.eta]
	return(list(mu=mu.vi, Gamma=Gamma.vi))
}
