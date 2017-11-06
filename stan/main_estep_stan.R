############################### Simulation - MCMC kernels (E-step) #############################
estep_stan<-function(kiter, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM) {
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
	post <- list(matrix(nrow = opt$nbiter.mcmc,ncol = ncol(phiM)))
	for (i in 1:(nrow(phiM))) {
		post[[i]] <- matrix(nrow = opt$nbiter.mcmc,ncol = ncol(phiM) )
	}

	
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]

	phiMc<-phiM
	for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau
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
	}
	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
	
	# Second stage
	
	if(opt$nbiter.mcmc[2]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1
		for (u in 1:opt$nbiter.mcmc[2]) {
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
		}
		varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
	}
	
	if(opt$nbiter.mcmc[3]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-kiter%%(nb.etas-1)+2
		if(is.nan(nrs2)) nrs2<-1 # to deal with case nb.etas=1
		for (u in 1:opt$nbiter.mcmc[3]) {
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
					post[[i]][u,] <- etaM[i,]
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
		}
		varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
	}

	if(opt$nbiter.mcmc[4]>0) {
		print('stan')
			for (i in 1:(nrow(phiM))) {
				Nj <- length(Dargs$yobs[Dargs$IdM==i])
				height <- Dargs$yobs[Dargs$IdM==i]
				age <- Dargs$XM[Dargs$IdM==i,]
				earn_dat <- list(N = Nj , #specify number of observations as a scalar
				                    height = height, # data vector
				                    age = age # data vector (predictor) 
				                    )

				earn_code <- 'data {
				  // First we declare all of our variables in the data block
				  int<lower=0> N;// Number of observations
				  vector[N] age; //Identify our predictor as a vector
				  vector[N] height;  //Identify our outcome variable as a vector
				}
				parameters {
				  vector[2] beta; //Our betas are a vector of length 2 (intercept and slope)
				  real<lower=0> sigma; //error parameter
				}
				model {
				  //Priors
				  beta[1] ~ normal( 1 , .001); //intercept
				  beta[2] ~ normal( 0 , 0.1 ); //slope
				  sigma ~ uniform( 0 , 1 ); //error
				  height ~ normal(beta[1] + beta[2] * age, sigma);
				}'

				browser()

				fit1 <- stan(model_code = earn_code, data = earn_dat,
				             warmup = 50,
				             iter = 500, 
				             chains = 4)
				
				fit1_samples = extract(fit1)
				str(fit1_samples)
				# subset just the betas
				betas = fit1_samples[[1]]

				phiM[i,]<-betas[end(betas)[1],]
			}

		
	}
	

	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	psiM<-transphi(phiM,Dargs$transform.par)
	fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
	if(Dargs$error.model=="exponential")
		fpred<-log(cutoff(fpred))
	gpred<-error(fpred,varList$pres)
	DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)^2+log(gpred)
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM, post = post))
}
