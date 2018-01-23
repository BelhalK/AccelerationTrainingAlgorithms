############################### Simulation - MCMC kernels (E-step) #############################
estep_cat_incremental<-function(kiter, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM,saemixObject,l,ind_rand) {
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
	# Uargs$nchains = 1
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
	psiM<-transphi(phiM,Dargs$transform.par)
	fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)


	DYF[Uargs$ind.ioM] <- -log(fpred)
	U.y<-colSums(DYF)
	
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]

	phiMc<-phiM
	map_range <- saemix.options$map.range
	
	index <- matrix(data=0,nrow=Dargs$NM,ncol=2)
	index[,1]<-1:Dargs$NM
	index[,2]<-rep(1:Dargs$N,Uargs$nchains)
	delete <- setdiff(1:Dargs$NM,index[index[,2] %in% l[ind_rand],1])

if (!(kiter %in% map_range)){
	# print('not in map range')
	for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau

		etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
		phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
		psiMc<-transphi(phiMc,Dargs$transform.par)
		fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
		DYF[Uargs$ind.ioM] <- -log(fpred)
		Uc.y<-colSums(DYF)
		deltau<-Uc.y-U.y
		deltau[delete] = 1000000
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
				DYF[Uargs$ind.ioM] <- -log(fpred)
				Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				# deltu[ind_rand] = 1000000
				deltu[delete] = 1000000
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
				DYF[Uargs$ind.ioM] <- -log(fpred)
				Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				# deltu[ind_rand] = 1000000
				deltu[delete] = 1000000
				ind<-which(deltu<(-log(runif(Dargs$NM))))
				etaM[ind,]<-etaMc[ind,]

				# for (i in 1:(nrow(phiM))) {
				# 	post[[i]][u,] <- etaM[i,]
				# }

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
}
	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))

	###############################################################################################
############   NEW KERNEl														############
###############################################################################################
if(opt$nbiter.mcmc[4]>0 & kiter %in% map_range) {
	# print('in map range')
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1

		#MAP calculation
		saemix.options<-saemixObject["options"]
	  	saemix.model<-saemixObject["model"]
	  	saemix.data<-saemixObject["data"]
	  	saemix.options$map <- TRUE
	  	saemixObject["results"]["omega"] <- omega.eta
	  	saemixObject["results"]["mean.phi"] <- mean.phi
	  	saemixObject["results"]["phi"] <- phiM
	  	# saemixObject["results"]["respar"] <- varList$pres

	  	i1.omega2<-saemixObject["model"]["indx.omega"]
	  	
	    iomega.phi1<-solve(saemixObject["results"]["omega"])
	  	id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
	  	xind<-saemixObject["data"]["data"][,saemixObject["data"]["name.predictors"], drop=FALSE]
	  	yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
	  	id.list<-unique(id)
	  	phi.map<-saemixObject["results"]["phi"]

	  	if (kiter %in% map_range){
	  	# print('start')
	  	for(i in 1:Dargs$NM) {
		    isuj<-id.list[i]
		    xi<-xind[id==isuj,,drop=FALSE]
		#    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
		    yi<-yobs[id==isuj]
		    idi<-rep(1,length(yi))
		    mean.phi1<-saemixObject["results"]["mean.phi"][i,i1.omega2]
		    phii<-saemixObject["results"]["phi"][i,]
		    phi1<-phii[i1.omega2]
		    
		  #      cond<-function(phi1){
				# 	conditional.distribution_cat(phi1,phii,idi,xi,yi,mean.phi1,i1.omega2,iomega.phi1,saemixObject["model"]["transform.par"],saemixObject["model"]["model"])
				# }
				
			 #    cma <- cmaNew()
				# cmaInit(cma,seed=395246,dimension=nb.etas,initialX=phii[i1.omega2])
				# res1 = cmaOptimDP(cma,cond,iterPrint=100);
				# mapi <- res1$xMat[res1$nIter,]
			 #    phi.map[i,i1.omega2]<-mapi
		    phi1.opti<-optim(par=phi1, fn=conditional.distribution_cat, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"])
		    # phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"],control = list(maxit = 2))
		    phi.map[i,i1.omega2]<-phi1.opti$par
		  }
		  # print('stop')
		 }
	  	map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
		map.psi<-data.frame(id=id.list,map.psi)
		map.phi<-data.frame(id=id.list,phi.map)

		psi_map <- as.matrix(map.psi[,-c(1)])
		phi_map <- as.matrix(map.phi[,-c(1)])
		eta_map <- phi_map[,varList$ind.eta] - mean.phiM[,varList$ind.eta]
		# eta_map <- phi_map - mean.phiM

		#gradient at the map estimation
		# gradf <- matrix(0L, nrow = length(fpred), ncol = nb.etas)
		gradp <- matrix(0L, nrow = Dargs$NM, ncol = nb.etas) 



		for (j in 1:nb.etas) {
			phi_map2 <- phi_map
			phi_map2[,j] <- phi_map[,j]+phi_map[,j]/100;
			psi_map2 <- transphi(phi_map2,saemixObject["model"]["transform.par"]) 
			fpred1<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
			DYF[Uargs$ind.ioM] <- -log(fpred1)
			l1 <- colSums(DYF)
			fpred2<-structural.model(psi_map2, Dargs$IdM, Dargs$XM)
			DYF[Uargs$ind.ioM] <- -log(fpred2)
			l2 <- colSums(DYF)
			for (i in 1:(Dargs$NM)){
				# r = 1:sum(Dargs$IdM == i)
    #             r = r+sum(as.matrix(gradf[,j]) != 0L)
				gradp[i,j] <- (l2[i] - l1[i])/(phi_map[i,j]/100)
			}
		}

		#calculation of the covariance matrix of the proposal
		
		# denom <- DYF
		fpred<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
		DYF[Uargs$ind.ioM] <- -log(fpred)
		denom <- colSums(DYF)

		
		Gamma <- chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
		z <- matrix(0L, nrow = length(fpred), ncol = 1) 
		for (i in 1:(Dargs$NM)){
			r = 1:sum(Dargs$IdM == i)
			r <- r+sum(as.matrix(z) != 0L)
            z[r] <- fpred[r]

			Gamma[[i]] <- solve(gradp[i,]%*%t(gradp[i,])/denom[i]^2+solve(omega.eta))
			# Gamma[[i]] <- omega.eta
			chol.Gamma[[i]] <- chol(Gamma[[i]])
			inv.Gamma[[i]] <- solve(Gamma[[i]])
		}
		

		etaM <- eta_map
		for (u in 1:opt$nbiter.mcmc[4]) {

			for(vk2 in 1:nb.etas) {
				etaMc<-etaM
				propc <- U.eta
				prop <- U.eta
				#generate candidate eta
					
				for (i in 1:(Dargs$NM)){
					Mi <- rnorm(nb.etas)%*%chol.Gamma[[i]]
					etaMc[i,]<- eta_map[i,] +Mi
				}


				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
				psiMc<-transphi(phiMc,Dargs$transform.par)
				fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
				DYF[Uargs$ind.ioM] <- -log(fpred)
				Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				

				for (i in 1:(Dargs$NM)){
					propc[i] <- 0.5*rowSums((etaMc[i,]-eta_map[i,])*(etaMc[i,]-eta_map[i,])%*%inv.Gamma[[i]])
					prop[i] <- 0.5*rowSums((etaM[i,]-eta_map[i,])*(etaM[i,]-eta_map[i,])%*%inv.Gamma[[i]])
				}



				deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind,]<-etaMc[ind,]
				
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM

			}
		}
	}

	post <- post_new <- 0

	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM, post = post, post_new=post_new))
}
