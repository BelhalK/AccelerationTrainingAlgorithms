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
	# map_range <- c(1:3)
	map_range <- saemix.options$map.range
	# map_range <- c(15,25,35)

	nb_replacement = round(saemix.options$nb.replacement*Dargs$NM/100)
	ind_rand = sample(1:Dargs$NM,(Dargs$NM-nb_replacement))

	if (!(kiter %in% map_range)){
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
		deltau[ind_rand] = 1000000
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
				deltu[ind_rand] = 1000000
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
				deltu[ind_rand] = 1000000
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
}
	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
	#New kernel with approx conditional

		if(opt$nbiter.mcmc[4]>0 & kiter %in% map_range) {
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
	  	saemixObject["results"]["respar"] <- varList$pres

	  	i1.omega2<-saemixObject["model"]["indx.omega"]
	    iomega.phi1<-solve(saemixObject["results"]["omega"][i1.omega2,i1.omega2])
	  	id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
	  	xind<-saemixObject["data"]["data"][,saemixObject["data"]["name.predictors"], drop=FALSE]
	  	yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
	  	id.list<-unique(id)
	  	phi.map<-saemixObject["results"]["phi"]


		K1 <- saemix.options$nbiter.saemix[1]
		K2 <- saemix.options$nbiter.saemix[2]
		K3 <- saemix.options$nbiter.saemix[3]

		K_gd <- 1
		gd_step <- saemix.options$step.gd
		if (kiter %in% map_range){

			if(kiter == 1){
				for(i in 1:Dargs$NM) {
			    isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			#    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-saemixObject["results"]["mean.phi"][i,i1.omega2]
			    phii<-saemixObject["results"]["phi"][i,]
			    phi1<-phii[i1.omega2]
			    phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"])
			    phi.map[i,i1.omega2]<-phi1.opti$par
			  }

			} else {


				for(i in ((1:Dargs$NM)-ind_rand)) {
			    isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			#    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-saemixObject["results"]["mean.phi"][i,i1.omega2]
			    phii<-saemixObject["results"]["phi"][i,]
			    phi1<-phii[i1.omega2]
			    phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"])
			    phi.map[i,i1.omega2]<-phi1.opti$par

			}
		  	
			 }
		else{
			map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
			map.psi<-data.frame(id=id.list,map.psi)
			map.phi<-data.frame(id=id.list,phi.map)

			psi_map <- as.matrix(map.psi[,-c(1)])
			phi_map <- as.matrix(map.phi[,-c(1)])
			eta_map <- phi_map - mean.phiM
			

			fpred1<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
			if(Dargs$error.model=="exponential")
				fpred1<-log(cutoff(fpred1))
			gpred1<-error(fpred1,varList$pres)
			DYF[Uargs$ind.ioM]<-1/sqrt(2*pi*gpred1)*exp(-0.5*((Dargs$yM-fpred1)/gpred1)**2)
			P.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
			P.eta<-0.5*rowSums(exp(-eta_map*(eta_map%*%somega)))


			for (k in 1:K_gd){
				gradp <- matrix(0L, nrow = Dargs$NM, ncol = nb.etas) 
			for(kj in 1:nb.etas){
				phi_map2 <- phi_map
				phi_map2[,kj] <- phi_map[,kj]+phi_map[,kj]/100
				psi_map2 <- transphi(phi_map2,Dargs$transform.par)
				eta_map2 <- phi_map2 - mean.phiM
				fpred2<-structural.model(psi_map2, Dargs$IdM, Dargs$XM)
				if(Dargs$error.model=="exponential")
					fpred2<-log(cutoff(fpred2))
				gpred2<-error(fpred2,varList$pres)
				DYF[Uargs$ind.ioM]<-1/sqrt(2*pi*gpred2)*exp(-0.5*((Dargs$yM-fpred2)/gpred2)**2)
				P2.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				P2.eta<-0.5*rowSums(exp(-eta_map2*(eta_map2%*%somega)))

				for (i in 1:Dargs$NM){
					gradp[i,kj] <- (P2.y[i]*P2.eta[i]-P.y[i]*P.eta[i])/(phi_map[i,kj]/100)
				}
			}
			phi.map <- phi.map + gd_step*gradp
			}
			
		}
	  	map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
		map.psi<-data.frame(id=id.list,map.psi)
		map.phi<-data.frame(id=id.list,phi.map)

		psi_map <- as.matrix(map.psi[,-c(1)])
		phi_map <- as.matrix(map.phi[,-c(1)])
		eta_map <- phi_map - mean.phiM
		
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
		

		#calculation of the covariance matrix of the proposal
	
		Gamma <- list(omega.eta,omega.eta)
		z <- matrix(0L, nrow = length(fpred), ncol = 1) 
		for (i in 1:(Dargs$NM)){
			r = 1:sum(Dargs$IdM == i)
			r <- r+sum(as.matrix(z) != 0L)
            z[r] <- gradf[r,1]
			Gamma[[i]] <- solve(t(gradf[r,])%*%gradf[r,]/(varList$pres[1])^2+solve(omega.eta))
		}
		# Gamma <- solve(t(gradf)%*%gradf/(varList$pres[1])^2+solve(omega.eta))
		# sGamma <- solve(Gamma)
		

		}
		for (u in 1:opt$nbiter.mcmc[4]) {
			for(vk2 in 1:nb.etas) {
				etaMc<-etaM
				etaM <- eta_map
				propc <- U.eta
				prop <- U.eta
				#generate candidate eta
					
				for (i in 1:(Dargs$NM)){
					M <- matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)%*%chol(Gamma[[i]])
					etaMc[i,]<- eta_map[i,] +M[i,]
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


				for (i in 1:(Dargs$NM)){
					propc[i] <- 0.5*rowSums((etaMc[i,]-eta_map[i,])*(etaMc[i,]-eta_map[i,])%*%solve(Gamma[[i]]))
					prop[i] <- 0.5*rowSums((etaM[i,]-eta_map[i,])*(etaM[i,]-eta_map[i,])%*%solve(Gamma[[i]]))
				}


				deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
				deltu[ind_rand] = 1000000
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind,]<-etaMc[ind,]
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


#New kernel with approx that hessian is null
	phiMold<-phiM
	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM[,]
	phiM[ind_rand,varList$ind.eta] <- phiMold[ind_rand,varList$ind.eta]

	# phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM, post = post))
}
