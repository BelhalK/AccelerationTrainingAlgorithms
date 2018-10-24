############################### Simulation - MCMC kernels (E-step) #############################

estep.incremental<-function(kiter, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM,saemixObject,l,ind_rand) {
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
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]

	ll <- sort(unique(l[ind_rand]),decreasing=FALSE)
	blockincremental <- NULL
	for (m in 1:Uargs$nchains){	
		blockincremental <- list.append(blockincremental,setdiff(1:Dargs$N, ll)+(m-1)*Dargs$N)
	}
	# chosen <- NULL
	# for (m in 1:Uargs$nchains){	
	# 	chosen <- list.append(chosen, ll+(m-1)*Dargs$N)
	# }
	chosen <- ll
	
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	phiMc<-phiM
	phiM[,varList$ind.eta] <- etaM[,varList$ind.eta] + mean.phiM[,varList$ind.eta]

	if(Dargs$type=="structural"){
		ll.y <- compute.LLy_c(phiM,varList$pres,Uargs,Dargs,DYF,1:Dargs$N,varList$omega)
		U.y<-ll.y$U
	} else{
		U.y <- compute.LLy_d(phiM,Uargs,Dargs,DYF)
	}
	
	# if (kiter==56) browser()
	# phiM <- etaM + mean.phiM
	# psiM<-transphi(phiM,Dargs$transform.par)
 #    tempsiM <- cbind(unique(Dargs$IdM), psiM)
 #    colnames(tempsiM) <- c("id",colnames(varList$omega))
 #    computePredictions(data.frame(tempsiM))[[1]]


	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	phiMc<-phiM
	for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau
		etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
		phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
		if(Dargs$type=="structural"){
			lly_c<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF,chosen,varList$omega)
			Uc.y<-lly_c$U
			indices <- lly_c$indices
			block <- unique(Dargs$IdM[indices])
		} else {
			Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
		}
		deltau<-Uc.y-U.y
		if (Dargs$monolix == TRUE){
		    deltau[block] <- 100000
		    deltau[blockincremental] <- 100000
		}
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
				if(Dargs$type=="structural"){
					lly_c<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF,chosen,varList$omega)
					Uc.y<-lly_c$U
					indices<-lly_c$indices
					block <- unique(Dargs$IdM[indices])
				} else {
					Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
				}
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				if (Dargs$monolix == TRUE){
				    deltu[block] <- 100000
				    deltu[blockincremental] <- 100000
				}
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
				if(Dargs$type=="structural"){
					lly_c<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF,chosen,varList$omega)
					Uc.y<-lly_c$U
					indices<-lly_c$indices
					block <- unique(Dargs$IdM[indices])
				} else {
					Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
				}
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				if (Dargs$monolix == TRUE){
					deltu[block] <- 100000
					deltu[blockincremental] <- 100000
				}
				ind<-which(deltu<(-log(runif(Dargs$NM))))
				etaM[ind,]<-etaMc[ind,]
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
	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM))
}
