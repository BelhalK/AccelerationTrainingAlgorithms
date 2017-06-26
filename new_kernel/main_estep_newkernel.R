############################### Simulation - MCMC kernels (E-step) #############################
estep_newkernel<-function(kiter, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM) {
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
	
	# "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)
	
	VK<-rep(c(1:nb.etas),2)
	Uargs$nchains = 1
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
	psiM<-transphi(phiM,Dargs$transform.par)
	fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
	if(Dargs$error.model=="exponential")
		fpred<-log(cutoff(fpred))
	gpred<-error(fpred,varList$pres)
	DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)^2+log(gpred)
	U.y<-colSums(DYF)
	


	post_rwm <- list(as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2)))
	for (i in 1:(nrow(phiM))) {
		post_rwm[[i]] <- as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2))
		names(post_rwm[[i]])[1] <- "iteration" 
		names(post_rwm[[i]])[ncol(post_rwm[[i]])] <- "individual"
		post_rwm[[i]][,1] <- 1:max(opt$nbiter.mcmc)
		post_rwm[[i]][,ncol(post_rwm[[i]])] <- i
	}


	

	post_vb <- list(as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2)))

	for (i in 1:(nrow(phiM))) {
		post_vb[[i]] <- as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2))
		names(post_vb[[i]])[1] <- "iteration" 
		names(post_vb[[i]])[ncol(post_vb[[i]])] <- "individual"
		post_vb[[i]][,1] <- 1:max(opt$nbiter.mcmc)
		post_vb[[i]][,ncol(post_vb[[i]])] <- i
	}

	post_vb_linear <- list(as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2)))

	for (i in 1:(nrow(phiM))) {
		post_vb_linear[[i]] <- as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2))
		names(post_vb_linear[[i]])[1] <- "iteration" 
		names(post_vb_linear[[i]])[ncol(post_vb_linear[[i]])] <- "individual"
		post_vb_linear[[i]][,1] <- 1:max(opt$nbiter.mcmc)
		post_vb_linear[[i]][,ncol(post_vb_linear[[i]])] <- i
	}

	
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]

	phiMc<-phiM

	for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau
		print(u)
		etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
		phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
		psiMc<-transphi(phiMc,Dargs$transform.par)
		fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
		if(Dargs$error.model=="exponential")
			fpred<-log(cutoff(fpred))
		gpred<-error(fpred,varList$pres)
		DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)^2+log(gpred)
		Uc.y<-colSums(DYF)
		deltau<-Uc.y-U.y
		ind<-which(deltau<(-1)*log(runif(Dargs$NM)))
		etaM[ind,]<-etaMc[ind,]
		U.y[ind]<-Uc.y[ind]

		U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1
		
			for(vk2 in 1:nb.etas) {
				etaMc<-etaM
				#				cat('vk2=',vk2,' nrs2=',nrs2,"\n")
				etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(Dargs$NM*nrs2), ncol=nrs2)%*%mydiag(varList$domega2[vk2,nrs2],nrow=1) # 2e noyau ? ou 1er noyau+permutation?
				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
				psiMc<-transphi(phiMc,Dargs$transform.par)
				fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred<-log(cutoff(fpred))
				gpred<-error(fpred,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
				Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind,]<-etaMc[ind,]
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM
			}
		
		varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))

		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-kiter%%(nb.etas-1)+2
		if(is.nan(nrs2)) nrs2<-1 # to deal with case nb.etas=1
		
			if(nrs2<nb.etas) {
				vk<-c(0,sample(c(1:(nb.etas-1)),nrs2-1))
				nb.iter2<-nb.etas
			} else {
				vk<-0:(nb.etas-1)
				#        if(nb.etas==1) vk<-c(0)
				nb.iter2<-1
			}
			for(k2 in 1:nb.iter2) {
				vk2<-VK[k2+vk]
				etaMc<-etaM
				etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(Dargs$NM*nrs2), ncol=nrs2)%*%mydiag(varList$domega2[vk2,nrs2])
				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
				psiMc<-transphi(phiMc,Dargs$transform.par)
				fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred<-log(cutoff(fpred))
				gpred<-error(fpred,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
				Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				ind<-which(deltu<(-log(runif(Dargs$NM))))
				etaM[ind,]<-etaMc[ind,]

				for (i in 1:(nrow(phiM))) {
					post_rwm[[i]][u,2:(ncol(post_rwm[[i]]) - 1)] <- etaM[i,]
				}


				#        if(kiter<20 | (kiter>150 & kiter<170)) {
				#        	cat("kiter=",kiter,length(ind),"  varList$ind.eta=",varList$ind.eta,"  nrs2=",nrs2,"\n")
				#        	print(head(etaMc))
				#        }
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM
			}

		varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))

	}


		#New kernel
		if(opt$nbiter.mcmc[4]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1

		#MAP calculation
		if(saemix.options$map) saemixObject<-map.saemix(saemixObject)
		browser()

		#gradient at the map estimation


		#calculation of the covariance matrix of the proposal
		Gamma <- solve(t(A)%*%A/(varList$pres[1])^2+solve(omega.eta))
		sGamma <- solve(Gamma)
		
		
		for (u in 1:opt$nbiter.mcmc[4]) {
			print(u)
			for(vk2 in 1:nb.etas) {
				etaMc<-etaM

				#generate candidate eta
				etaMc<- eta_map +matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)%*%chol(Gamma)

				
				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
				psiMc<-transphi(phiMc,Dargs$transform.par)
				fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred<-log(cutoff(fpred))
				gpred<-error(fpred,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
				Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))

				propc <- 0.5*rowSums((etaMc-eta_map)*(etaMc-eta_map)%*%sGamma))
				prop <- 0.5*rowSums((etaM-eta_map)*(etaM-eta_map)%*%sGamma))
				
				deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind,]<-etaMc[ind,]
				for (i in 1:(nrow(phiM))) {
					post_vb_linear[[i]][u,2:(ncol(post_vb_linear[[i]]) - 1)] <- etaM[i,]
				}
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM


				# #Or Use the output of VI as the posterior distrib we simulate from
				# etaM[ind,]<-etaMc[ind,]
				# for (i in 1:(nrow(phiM))) {
				# 	post_vb[[i]][u,2:(ncol(post_vb[[i]]) - 1)] <- etaM[i,]
				# }

			}
		}
	}

	
	
		
	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM, post_rwm = post_rwm,post_vb = post_vb,post_vb_linear = post_vb_linear))
}
