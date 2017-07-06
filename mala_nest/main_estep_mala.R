############################### Simulation - MCMC kernels (E-step) #############################
estep_mala<-function(kiter, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM,saemixObject) {
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

	post_mala <- list(as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2)))

	for (i in 1:(nrow(phiM))) {
		post_mala[[i]] <- as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2))
		names(post_mala[[i]])[1] <- "iteration" 
		names(post_mala[[i]])[ncol(post_mala[[i]])] <- "individual"
		post_mala[[i]][,1] <- 1:max(opt$nbiter.mcmc)
		post_mala[[i]][,ncol(post_mala[[i]])] <- i
	}


	x <- list(as.data.frame(matrix(nrow = Dargs$NM,ncol = ncol(phiM))))

	for (i in 1:(nrow(phiM))) {
		x[[i]] <- as.data.frame(matrix(nrow = Dargs$NM,ncol = ncol(phiM)))
	}


	
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]

	phiMc<-phiM

	for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau
		# print(u)
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


		#MALA
		if(opt$nbiter.mcmc[4]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1
		adap <- rep(1, Dargs$NM)
		sigma <- 0.01
		gamma <- 0.01
		
		for (u in 1:opt$nbiter.mcmc[4]) {
			# print(u)
			
			etaMc<-etaM
			propc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			prop <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradU <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradUc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			#Gradient in current eta
			phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
			psiM<-transphi(phiM,Dargs$transform.par)
			fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred<-log(cutoff(fpred))
			gpred<-error(fpred,varList$pres)
			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
			U.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
			U.eta<-0.5*rowSums(etaM*(etaM%*%somega))

			for (kj in 1:(nb.etas)){
				etaM2 <- etaM
				phiM2 <- phiM
				etaM2[,kj] <- etaM[,kj] + etaM[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred2<-log(cutoff(fpred2))
				gpred2<-error(fpred2,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred2)/gpred2)**2+log(gpred2)
				U2.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				
				for (i in 1:Dargs$NM){
					gradU[i,kj] <- -(U2.y[i]-U.y[i]+U2.eta[i]-U.eta[i])/(etaM[i,kj]/100)
				}
			}
			# 

			if (u>1){
				adap <- adap - gamma*(deltu + log(0.57))
			}
			
			Z <- matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)

			for (i in 1:Dargs$NM){
				etaMc[i,] <- etaM[i,] + sigma*adap[i]*gradU[i,] + sqrt(2*sigma*adap[i])*Z[i,]
			}
			

			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
			psiMc<-transphi(phiMc,Dargs$transform.par)
			fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred<-log(cutoff(fpred))
			gpred<-error(fpred,varList$pres)
			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
			Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
			Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))

			#Gradient in candidate eta

			for (kj in 1:(nb.etas)){
				etaM2 <- etaMc
				phiM2 <- phiMc
				etaM2[,kj] <- etaMc[,kj] + etaMc[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred2<-log(cutoff(fpred2))
				gpred2<-error(fpred2,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred2)/gpred2)**2+log(gpred2)
				U2.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				for (i in 1:Dargs$NM){
					gradUc[i,kj] <- -(U2.y[i]-Uc.y[i]+U2.eta[i]-Uc.eta[i])/(etaMc[i,kj]/100)
				}
			}

			
			for (i in 1:(Dargs$NM)){
				propc[i,] <- ((etaMc[i,]-etaM[i,] - sigma*adap[i]*gradU[i,])/sqrt(2*sigma*adap[i]))^2
				prop[i,] <- ((etaM[i,]-etaMc[i,] - sigma*adap[i]*gradUc[i,])/sqrt(2*sigma*adap[i]))^2
			}
			

			P<-0.5*rowSums(prop)
			Pc<-0.5*rowSums(propc)

			deltu<-Uc.y-U.y+Uc.eta-U.eta + P - Pc
			ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
			# print(length(ind)/Dargs$NM)
			etaM[ind,]<-etaMc[ind,]
			for (i in 1:(nrow(phiM))) {
				post_mala[[i]][u,2:(ncol(post_mala[[i]]) - 1)] <- etaM[i,]
			}
			U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
			U.eta[ind]<-Uc.eta[ind]
			nbc2<-nbc2+length(ind)
			nt2<-nt2+Dargs$NM


			# #Or Use the output of VI as the posterior distrib we simulate from
			# etaM[ind,]<-etaMc[ind,]
			# for (i in 1:(nrow(phiM))) {
			# 	post_vb[[i]][u,2:(ncol(post_vb[[i]]) - 1)] <- etaM[i,]
			# }

			
		}
	}

	#Nesterov
		if(opt$nbiter.mcmc[5]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1
		adap <- rep(1, Dargs$NM)
		sigma <- 0.01
		gamma <- 0.01
		
		for (u in 1:opt$nbiter.mcmc[5]) {
			# print(u)
			
			etaMc<-etaM
			propc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			prop <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradU <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradUc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			#Gradient in current eta
			phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
			psiM<-transphi(phiM,Dargs$transform.par)
			fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred<-log(cutoff(fpred))
			gpred<-error(fpred,varList$pres)
			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
			U.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
			U.eta<-0.5*rowSums(etaM*(etaM%*%somega))

			for (kj in 1:(nb.etas)){
				etaM2 <- etaM
				phiM2 <- phiM
				etaM2[,kj] <- etaM[,kj] + etaM[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred2<-log(cutoff(fpred2))
				gpred2<-error(fpred2,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred2)/gpred2)**2+log(gpred2)
				U2.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				
				for (i in 1:Dargs$NM){
					gradU[i,kj] <- -(U2.y[i]-U.y[i]+U2.eta[i]-U.eta[i])/(etaM[i,kj]/100)
				}
			}
			# 

			if (u>1){
				adap <- adap - gamma*(deltu + log(0.57))
			}
			
			Z <- matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)

			R=0.55
			a<-1
			if (u>2){

				if (u<100){
					a <- 0.5
					for (i in 1:Dargs$NM){
						etaMc[i,] <- etaM[i,] + sigma*adap[i]*gradU[i,] +R*(etaM[i,] - x[[u-2]][i,]) + sqrt(2*a*sigma*adap[i])*Z[i,]
					}
					a<-1
				}
				else {
					a <- 1
					for (i in 1:Dargs$NM){
						etaMc[i,] <- etaM[i,] + sigma*adap[i]*gradU[i,] +R*(etaM[i,] - x[[u-2]][i,]) + sqrt(2*sigma*adap[i])*Z[i,]
					}
				}

				
			}
			

			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
			psiMc<-transphi(phiMc,Dargs$transform.par)
			fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred<-log(cutoff(fpred))
			gpred<-error(fpred,varList$pres)
			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
			Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
			Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))

			#Gradient in candidate eta

			for (kj in 1:(nb.etas)){
				etaM2 <- etaMc
				phiM2 <- phiMc
				etaM2[,kj] <- etaMc[,kj] + etaMc[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred2<-log(cutoff(fpred2))
				gpred2<-error(fpred2,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred2)/gpred2)**2+log(gpred2)
				U2.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				for (i in 1:Dargs$NM){
					gradUc[i,kj] <- -(U2.y[i]-Uc.y[i]+U2.eta[i]-Uc.eta[i])/(etaMc[i,kj]/100)
				}
			}

			
			for (i in 1:(Dargs$NM)){
				propc[i,] <- ((etaMc[i,]-etaM[i,] - sigma*adap[i]*gradU[i,])/sqrt(2*a*sigma*adap[i]))^2
				prop[i,] <- ((etaM[i,]-etaMc[i,] - sigma*adap[i]*gradUc[i,])/sqrt(2*a*sigma*adap[i]))^2
			}
			

			P<-0.5*rowSums(prop)
			Pc<-0.5*rowSums(propc)

			deltu<-Uc.y-U.y+Uc.eta-U.eta + P - Pc
			ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
			# print(length(ind)/Dargs$NM)
			etaM[ind,]<-etaMc[ind,]
			x[[u]] <- etaM
			for (i in 1:(nrow(phiM))) {
				post_vb[[i]][u,2:(ncol(post_vb[[i]]) - 1)] <- etaM[i,]
			}
			U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
			U.eta[ind]<-Uc.eta[ind]
			nbc2<-nbc2+length(ind)
			nt2<-nt2+Dargs$NM

			
		}
	}

	#AMALA
		if(opt$nbiter.mcmc[6]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1
		adap <- rep(1, Dargs$NM)
		sigma <- 0.01
		gamma <- 0.01
		
		for (u in 1:opt$nbiter.mcmc[6]) {
			# print(u)
			
			etaMc<-etaM
			propc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			prop <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradU <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradUc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			#Gradient in current eta
			phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
			psiM<-transphi(phiM,Dargs$transform.par)
			fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred<-log(cutoff(fpred))
			gpred<-error(fpred,varList$pres)
			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
			U.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
			U.eta<-0.5*rowSums(etaM*(etaM%*%somega))

			for (kj in 1:(nb.etas)){
				etaM2 <- etaM
				phiM2 <- phiM
				etaM2[,kj] <- etaM[,kj] + etaM[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred2<-log(cutoff(fpred2))
				gpred2<-error(fpred2,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred2)/gpred2)**2+log(gpred2)
				U2.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				
				for (i in 1:Dargs$NM){
					gradU[i,kj] <- -(U2.y[i]-U.y[i]+U2.eta[i]-U.eta[i])/(etaM[i,kj]/100)
				}
			}
			# 

			if (u>1){
				adap <- adap - gamma*(deltu + log(0.57))
			}
			

			Z <- matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)

			sd <- sigma/max(sigma,norm(gradU))*gradU
			G <- t(sd)%*%sd
			G <- G+10*diag(nb.etas)
			sG <- solve(G)
			sd <- chol(G)

			B <- Z%*%t(sd)
			# B <- sigma/max(sigma,norm(gradU))*gradU*Z
			
			
			
			for (i in 1:Dargs$NM){
				# etaMc[i,] <- etaM[i,] + sigma/max(sigma,norm(gradU))*gradU[i,] + sqrt(2*sigma*adap[i])*Z[i,]
				etaMc[i,] <- etaM[i,] + sigma*sigma/max(sigma,norm(gradU))*gradU[i,] + sqrt(sigma)*B[i,]
			}
			

			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
			psiMc<-transphi(phiMc,Dargs$transform.par)
			fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred<-log(cutoff(fpred))
			gpred<-error(fpred,varList$pres)
			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
			Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
			Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))

			#Gradient in candidate eta

			for (kj in 1:(nb.etas)){
				etaM2 <- etaMc
				phiM2 <- phiMc
				etaM2[,kj] <- etaMc[,kj] + etaMc[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred2<-log(cutoff(fpred2))
				gpred2<-error(fpred2,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred2)/gpred2)**2+log(gpred2)
				U2.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				for (i in 1:Dargs$NM){
					gradUc[i,kj] <- -(U2.y[i]-Uc.y[i]+U2.eta[i]-Uc.eta[i])/(etaMc[i,kj]/100)
				}
			}

			
			# browser()

			for (i in 1:(Dargs$NM)){
				# propc[i,] <- ((etaMc[i,]-etaM[i,] - sigma*adap[i]*gradU[i,])/sqrt(2*sigma*adap[i]))^2
				# prop[i,] <- ((etaM[i,]-etaMc[i,] - sigma*adap[i]*gradUc[i,])/sqrt(2*sigma*adap[i]))^2

				propc[i,] <- ((etaMc[i,]-etaM[i,] - sigma/max(sigma,norm(gradU))*gradU[i,])/sqrt(sigma))^2%*%sG
				prop[i,] <- ((etaM[i,]-etaMc[i,] - sigma/max(sigma,norm(gradU))*gradUc[i,])/sqrt(sigma))^2%*%sG
			}
			

			P<-0.5*rowSums(prop)
			Pc<-0.5*rowSums(propc)

			deltu<-Uc.y-U.y+Uc.eta-U.eta + P - Pc
			ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
			# print(length(ind)/Dargs$NM)
			etaM[ind,]<-etaMc[ind,]
			for (i in 1:(nrow(phiM))) {
				post_mala[[i]][u,2:(ncol(post_mala[[i]]) - 1)] <- etaM[i,]
			}
			U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
			U.eta[ind]<-Uc.eta[ind]
			nbc2<-nbc2+length(ind)
			nt2<-nt2+Dargs$NM


			# #Or Use the output of VI as the posterior distrib we simulate from
			# etaM[ind,]<-etaMc[ind,]
			# for (i in 1:(nrow(phiM))) {
			# 	post_vb[[i]][u,2:(ncol(post_vb[[i]]) - 1)] <- etaM[i,]
			# }

			
		}
	}	
	
	#Non Reversible
		if(opt$nbiter.mcmc[7]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1
		adap <- rep(1, Dargs$NM)
		sigma <- 0.01
		gamma <- 0.01
		
		for (u in 1:opt$nbiter.mcmc[7]) {
			# print(u)
			
			etaMc<-etaM
			propc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			prop <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradU <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradUc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			#Gradient in current eta
			phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
			psiM<-transphi(phiM,Dargs$transform.par)
			fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred<-log(cutoff(fpred))
			gpred<-error(fpred,varList$pres)
			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
			U.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
			U.eta<-0.5*rowSums(etaM*(etaM%*%somega))

			for (kj in 1:(nb.etas)){
				etaM2 <- etaM
				phiM2 <- phiM
				etaM2[,kj] <- etaM[,kj] + etaM[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred2<-log(cutoff(fpred2))
				gpred2<-error(fpred2,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred2)/gpred2)**2+log(gpred2)
				U2.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				
				for (i in 1:Dargs$NM){
					gradU[i,kj] <- -(U2.y[i]-U.y[i]+U2.eta[i]-U.eta[i])/(etaM[i,kj]/100)
				}
			}
			# 

			if (u>1){
				adap <- adap - gamma*(deltu + log(0.57))
			}
			
			Z <- matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)

			for (i in 1:Dargs$NM){
				etaMc[i,] <- etaM[i,] + sigma*adap[i]*gradU[i,] + sqrt(2*sigma*adap[i])*Z[i,]
			}
			

			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
			psiMc<-transphi(phiMc,Dargs$transform.par)
			fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred<-log(cutoff(fpred))
			gpred<-error(fpred,varList$pres)
			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
			Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
			Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))

			#Gradient in candidate eta

			for (kj in 1:(nb.etas)){
				etaM2 <- etaMc
				phiM2 <- phiMc
				etaM2[,kj] <- etaMc[,kj] + etaMc[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred2<-log(cutoff(fpred2))
				gpred2<-error(fpred2,varList$pres)
				DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred2)/gpred2)**2+log(gpred2)
				U2.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				for (i in 1:Dargs$NM){
					gradUc[i,kj] <- -(U2.y[i]-Uc.y[i]+U2.eta[i]-Uc.eta[i])/(etaMc[i,kj]/100)
				}
			}

			
			for (i in 1:(Dargs$NM)){
				propc[i,] <- ((etaMc[i,]-etaM[i,] - sigma*adap[i]*gradU[i,])/sqrt(2*sigma*adap[i]))^2
				prop[i,] <- ((etaM[i,]-etaMc[i,] - sigma*adap[i]*gradUc[i,])/sqrt(2*sigma*adap[i]))^2
			}
			

			P<-0.5*rowSums(prop)
			Pc<-0.5*rowSums(propc)

			deltu<-Uc.y-U.y+Uc.eta-U.eta + P - Pc
			ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
			# print(length(ind)/Dargs$NM)
			etaM[ind,]<-etaMc[ind,]
			for (i in 1:(nrow(phiM))) {
				post_mala[[i]][u,2:(ncol(post_mala[[i]]) - 1)] <- etaM[i,]
			}
			U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
			U.eta[ind]<-Uc.eta[ind]
			nbc2<-nbc2+length(ind)
			nt2<-nt2+Dargs$NM


			# #Or Use the output of VI as the posterior distrib we simulate from
			# etaM[ind,]<-etaMc[ind,]
			# for (i in 1:(nrow(phiM))) {
			# 	post_vb[[i]][u,2:(ncol(post_vb[[i]]) - 1)] <- etaM[i,]
			# }

			
		}
	}	

	
	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM, post_rwm = post_rwm,post_vb = post_vb,post_mala = post_mala))
}

