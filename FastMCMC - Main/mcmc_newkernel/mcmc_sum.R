mcmc_sum<-function(model,data,control=list(),iter) {
	# E-step - simulate unknown parameters
	# Input: kiter, Uargs, structural.model, mean.phi (unchanged)
	# Output: varList, DYF, phiM (changed)
	



# create progress bar
pb <- txtProgressBar(min = 0, max = iter, style = 1)


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
	Uargs$nchains = 1
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))


	#MAP calculation

	saemix.options<-saemixObject["options"]
  	saemix.model<-saemixObject["model"]
  	saemix.data<-saemixObject["data"]
  	saemix.options$map <- TRUE
  	saemixObject["results"]["omega"] <- omega.eta
  	saemixObject["results"]["mean.phi"] <- mean.phi
  	saemixObject["results"]["phi"] <- phiM
  	saemixObject["results"]["respar"] <- varList$pres

  	i1.omega2<-saemixObject["model"]["indx.omega"]
    iomega.phi1<-solve(saemixObject["results"]["omega"][i1.omega2,i1.omega2])
  	id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  	xind<-saemixObject["data"]["data"][,saemixObject["data"]["name.predictors"], drop=FALSE]
  	yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  	id.list<-unique(id)
  	phi.map<-saemixObject["results"]["phi"]

  	for(i in 1:Dargs$NM) {
	    isuj<-id.list[i]
	    xi<-xind[id==isuj,,drop=FALSE]
	#    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
	    yi<-yobs[id==isuj]
	    idi<-rep(1,length(yi))
	    mean.phi1<-mean.phi
	    phii<-phiM[i,]
	    phi1<-phii[i1.omega2]
	    phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"])
	    # phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"],control = list(maxit = 2))
	    phi.map[i,i1.omega2]<-phi1.opti$par
	  }

  	map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
	map.psi<-data.frame(id=id.list,map.psi)
	map.phi<-data.frame(id=id.list,phi.map)

	psi_map <- as.matrix(map.psi[,-c(1)])
	phi_map <- as.matrix(map.phi[,-c(1)])
	eta_map <- phi_map - mean.phiM



	etaM <- eta_map
	phiM <- phi_map
	psiM<-transphi(phiM,Dargs$transform.par)
	fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
	gpred<-error(fpred,varList$pres)
	DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)^2+log(gpred)
	U.y<-colSums(DYF)
	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
	phiMc<-phiM



	eta_list <- list(as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2)))

	for (i in 1:(nrow(phiM))) {
		eta_list[[i]] <- as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = ncol(phiM)+2))
		names(eta_list[[i]])[1] <- "iteration" 
		names(eta_list[[i]])[ncol(eta_list[[i]])] <- "individual"
		eta_list[[i]][,1] <- 1:max(opt$nbiter.mcmc)
		eta_list[[i]][,ncol(eta_list[[i]])] <- i
	}

	dens_Ueta <- list(as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = nb.etas)))
	for (i in 1:(nrow(phiM))) {
		dens_Ueta[[i]] <- as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = nb.etas))
		names(dens_Ueta[[i]])[1] <- "iteration" 
		names(dens_Ueta[[i]])[ncol(dens_Ueta[[i]])] <- "individual"
		dens_Ueta[[i]][,1] <- 1:max(opt$nbiter.mcmc)
		dens_Ueta[[i]][,ncol(dens_Ueta[[i]])] <- i
	}

	densy <- list(as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = nb.etas)))
	for (i in 1:(nrow(phiM))) {
		densy[[i]] <- as.data.frame(matrix(nrow = max(opt$nbiter.mcmc),ncol = nb.etas))
		names(densy[[i]])[1] <- "iteration" 
		names(densy[[i]])[ncol(densy[[i]])] <- "individual"
		densy[[i]][,1] <- 1:max(opt$nbiter.mcmc)
		densy[[i]][,ncol(densy[[i]])] <- i
	}

for(u in 1:opt$nbiter.mcmc[1]) {


		Sys.sleep(0.0001)
   			# update progress bar
			setTxtProgressBar(pb,u)

			for (i in 1:(nrow(phiM))) {
				eta_list[[i]][u,2:(nb.etas+2 - 1)] <- etaM[i,]
			}
			for (i in 1:(nrow(phiM))) {
				dens_Ueta[[i]][u,2] <- U.eta[i]
			}
			for (i in 1:(nrow(phiM))) {
				densy[[i]][u,2] <- U.y[i]
			}
		etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
		phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
		psiMc<-transphi(phiMc,Dargs$transform.par)
		fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
		if(Dargs$error.model=="exponential")
			fpred<-log(cutoff(fpred))
		gpred<-error(fpred,varList$pres)
		DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)^2+log(gpred)
		Uc.y<-colSums(DYF)
		Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
		deltau<-Uc.y-U.y
		ind<-which(deltau<(-1)*log(runif(Dargs$NM)))
		etaM[ind,]<-etaMc[ind,]
		U.y[ind]<-Uc.y[ind]
		U.eta[ind]<-Uc.eta[ind]

	if(opt$nbiter.mcmc[2]>0) {
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
}

if(opt$nbiter.mcmc[3]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1%%(nb.etas-1)+2
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
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]

				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM
			}
			

			varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
		}
	}



		#New kernel
	if(opt$nbiter.mcmc[4]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1

		
		#gradient at the map estimation
		gradf <- matrix(0L, nrow = length(fpred), ncol = nb.etas) 

		for (j in 1:nb.etas) {
			phi_map2 <- phi_map
			phi_map2[,j] <- phi_map[,j]+phi_map[,j]/100;
			psi_map2 <- transphi(phi_map2,saemixObject["model"]["transform.par"]) 
			fpred1<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
			fpred2<-structural.model(psi_map2, Dargs$IdM, Dargs$XM)
			for (i in 1:(Dargs$NM)){
				r = 1:sum(Dargs$IdM == i)
                r = r+sum(as.matrix(gradf[,j]) != 0L)
				gradf[r,j] <- (fpred2[r] - fpred1[r])/(phi_map[i,j]/100)
			}
		}

		#calculation of the covariance matrix of the proposal
	
		Gamma <- chol.Gamma <- inv.Gamma <- inv.sum <- list(omega.eta,omega.eta)
		z <- matrix(0L, nrow = length(fpred), ncol = 1) 
		for (i in 1:(Dargs$NM)){
			r = 1:sum(Dargs$IdM == i)
			r <- r+sum(as.matrix(z) != 0L)
            z[r] <- gradf[r,1]
			Gamma[[i]] <- solve(t(gradf[r,])%*%gradf[r,]/(varList$pres[1])^2+solve(omega.eta))
			chol.Gamma[[i]] <- chol(Gamma[[i]])
			inv.Gamma[[i]] <- solve(Gamma[[i]])
			inv.sum[[i]] <- solve(Gamma[[i]]/4+varList$domega2/4)
		}
		
		etaMc<-etaM
		etaM <- eta_map
		propc <- U.eta
		prop <- U.eta
		for (u in 1:opt$nbiter.mcmc[4]) {
			Sys.sleep(0.0001)
   			# update progress bar
			setTxtProgressBar(pb,u)

			for (i in 1:(nrow(phiM))) {
				eta_list[[i]][u,2:(nb.etas+2 - 1)] <- etaM[i,]
			}
			for (i in 1:(nrow(phiM))) {
				dens_Ueta[[i]][u,2] <- U.eta[i]
			}
			for (i in 1:(nrow(phiM))) {
				densy[[i]][u,2] <- U.y[i]
			}
				#generate candidate eta
				

				nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
				nrs2<-1
		
				for(vk2 in 1:nb.etas) {
					etaMc<-etaM
					#				cat('vk2=',vk2,' nrs2=',nrs2,"\n")
					for (i in 1:(Dargs$NM)){
					Mi <- rnorm(nb.etas)%*%chol.Gamma[[i]]
					etaMc[i,vk2]<-0.5*(etaM[i,vk2]+matrix(rnorm(Dargs$NM*nrs2), ncol=nrs2)%*%mydiag(varList$domega2[vk2,nrs2],nrow=1))+0.5*(eta_map[i,vk2] +Mi[,vk2])
					}
					
					phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
					psiMc<-transphi(phiMc,Dargs$transform.par)
					fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
					gpred<-error(fpred,varList$pres)
					DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
					Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
					Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))

					for (i in 1:(Dargs$NM)){
						propc[i] <- 0.5*rowSums((etaMc[i,]-eta_map[i,]/2-etaM[i,]/2)*(etaMc[i,]-eta_map[i,]/2-etaM[i,]/2)%*%inv.sum[[i]])
						prop[i] <- 0.5*rowSums((etaM[i,]-eta_map[i,]/2-etaM[i,]/2)*(etaM[i,]-eta_map[i,]/2-etaM[i,]/2)%*%inv.sum[[i]])
					}

					deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
					ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
					etaM[ind,]<-etaMc[ind,]
					U.y[ind]<-Uc.y[ind]
					U.eta[ind]<-Uc.eta[ind]
				
					nbc2[vk2]<-nbc2[vk2]+length(ind)
					nt2[vk2]<-nt2[vk2]+Dargs$NM
				}
		
		varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))


			}
		
	}
	
	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	close(pb)
	return(list(eta = eta_list,denseta = dens_Ueta,densy = densy))
}

