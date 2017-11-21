############################### Simulation - MCMC kernels (E-step) #############################
estep_time_mamyula<-function(kiter, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM,saemixObject) {
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
	DYF[Uargs$ind.ioM]<- -fpred
	U.y<-colSums(DYF)

	# U.y <- -fpred

	post <- list(matrix(nrow = opt$nbiter.mcmc,ncol = ncol(phiM)))
	for (i in 1:(nrow(phiM))) {
		post[[i]] <- matrix(nrow = opt$nbiter.mcmc,ncol = ncol(phiM) )
	}

	
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]

	phiMc<-phiM
	map_range <- saemix.options$map.range

if (!(kiter %in% map_range)){
	
	for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau

		etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
		phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
		psiMc<-transphi(phiMc,Dargs$transform.par)
		fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
		DYF[Uargs$ind.ioM]<- -fpred
		Uc.y<-colSums(DYF)
		# if (kiter>25){browser()}
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
				
				DYF[Uargs$ind.ioM]<- -fpred
				Uc.y<-colSums(DYF)
				# Uc.y <- -fpred
				# Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
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
				DYF[Uargs$ind.ioM]<- -fpred
				Uc.y<-colSums(DYF)
				# Uc.y <- -fpred
				
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
}
	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))

	###############################################################################################
############   NEW KERNEl														############
###############################################################################################
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
	  	# saemixObject["results"]["respar"] <- varList$pres

	  	i1.omega2<-saemixObject["model"]["indx.omega"]
	    iomega.phi1<-solve(saemixObject["results"]["omega"][i1.omega2,i1.omega2])
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
		    phi1.opti<-optim(par=phi1, fn=conditional.distribution_time, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"])
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
		
		#gradient at the map estimation
		# gradf <- matrix(0L, nrow = length(fpred), ncol = nb.etas)
		gradp <- matrix(0L, nrow = Dargs$NM, ncol = nb.etas) 

		# for (j in 1:nb.etas) {
		# 	phi_map2 <- phi_map
		# 	phi_map2[,j] <- phi_map[,j]+phi_map[,j]/100;
		# 	psi_map2 <- transphi(phi_map2,saemixObject["model"]["transform.par"]) 
		# 	fpred1<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
		# 	fpred2<-structural.model(psi_map2, Dargs$IdM, Dargs$XM)
		# 	for (i in 1:(Dargs$NM)){
		# 		r = 1:sum(Dargs$IdM == i)
  #               r = r+sum(as.matrix(gradf[,j]) != 0L)
		# 		gradf[r,j] <- (fpred2[r] - fpred1[r])/(phi_map[i,j]/100)
		# 	}
		# }

		for (j in 1:nb.etas) {
			phi_map2 <- phi_map
			phi_map2[,j] <- phi_map[,j]+phi_map[,j]/100;
			psi_map2 <- transphi(phi_map2,saemixObject["model"]["transform.par"]) 
			fpred1<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
			DYF[Uargs$ind.ioM]<- fpred1
			l1<-colSums(DYF)
			fpred2<-structural.model(psi_map2, Dargs$IdM, Dargs$XM)
			DYF[Uargs$ind.ioM]<- fpred2
			l2<-colSums(DYF)

			for (i in 1:(Dargs$NM)){
				# r = 1:sum(Dargs$IdM == i)
    #             r = r+sum(as.matrix(gradf[,j]) != 0L)
				gradp[i,j] <- (l2[i] - l1[i])/(phi_map[i,j]/100)
			}
		}

		#calculation of the covariance matrix of the proposal
		
		# denom <- DYF
		fpred<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
		
		DYF[Uargs$ind.ioM]<- fpred
		denom <- colSums(DYF)

		
		Gamma <- chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
		z <- matrix(0L, nrow = length(fpred), ncol = 1) 
		for (i in 1:(Dargs$NM)){
			# r = 1:sum(Dargs$IdM == i)
			# r <- r+sum(as.matrix(z) != 0L)
   #          z[r] <- fpred[r]

			Gamma[[i]] <- solve(gradp[i,]%*%t(gradp[i,])/denom[i]^2+solve(omega.eta))
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
				DYF[Uargs$ind.ioM]<- -fpred
				Uc.y<-colSums(DYF)
				# Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))


				for (i in 1:(Dargs$NM)){
					propc[i] <- 0.5*rowSums((etaMc[i,]-eta_map[i,])*(etaMc[i,]-eta_map[i,])%*%inv.Gamma[[i]])
					prop[i] <- 0.5*rowSums((etaM[i,]-eta_map[i,])*(etaM[i,]-eta_map[i,])%*%inv.Gamma[[i]])
				}


				deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind]<-etaMc[ind]
				# for (i in 1:(nrow(phiM))) {
				# 	post_newkernel[[i]][u,2:(ncol(post_newkernel[[i]]) - 1)] <- etaM[i]
				# }
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM


				# #Or Use the output of VI as the posterior distrib we simulate from
				# etaM[ind,]<-etaMc[ind,]
				# for (i in 1:(nrow(phiM))) {
				# 	post_vb[[i]][u,2:(ncol(post_vb[[i]]) - 1)] <- etaM[i]
				# }

			}
		}
	}

	if(opt$nbiter.mcmc[5]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1
		adap <- rep(1, Dargs$NM)
		sigma <- saemix.options$sigma.val
		gamma <- saemix.options$gamma.val

		acc <- 0
		lambda<-saemix.options$lambda.val
		for (u in 1:opt$nbiter.mcmc[5]) {
			# print(u)
			
##### find the prox of g (unique minimizer of a regularized g=p(y|z))
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
		  	phi.prox<-saemixObject["results"]["phi"]

		  	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
		  	
		  	
		  	for(i in 1:Dargs$NM) {
			    isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			#    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-mean.phiM
			    # phii<-saemixObject["results"]["phi"][i,]
			    phii<-phiM
			    phi1<-phii[i1.omega2]
			    
			    phi1.opti<-optim(par=phi1, fn=proximal.function, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1,
			     trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], 
			     pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"],
			     lambda=lambda,current_phi=phiM,,control = list(maxit =1))
			    # phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"],control = list(maxit = 2))
			    phi.prox[i,i1.omega2]<-phi1.opti$par
			  }

			prox.psi<-transphi(phi.prox,saemixObject["model"]["transform.par"])
			prox.psi<-data.frame(id=id.list,prox.psi)
			prox.phi<-data.frame(id=id.list,phi.prox)

			psi_prox <- as.matrix(prox.psi[,-c(1)])
			phi_prox <- as.matrix(prox.phi[,-c(1)])
			eta_prox <- phi_prox - mean.phiM


			etaMc<-etaM
			propc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			prop <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradU <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			gradUc <- matrix(nrow = Dargs$NM,ncol = nb.etas)
			#Gradient in current eta
			phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
			psiM<-transphi(phiM,Dargs$transform.par)
			fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
			DYF[Uargs$ind.ioM]<- -fpred
			U.y<-colSums(DYF)
			U.eta<-0.5*rowSums(etaM*(etaM%*%somega))

			for (kj in 1:(nb.etas)){
				etaM2 <- etaM
				phiM2 <- phiM
				etaM2[,kj] <- etaM[,kj] + etaM[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				DYF[Uargs$ind.ioM]<- -fpred2
				U2.y<-colSums(DYF)
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				
				# for (i in 1:Dargs$NM){
				# 	# gradU[i,kj] <- -(U2.y[i]-U.y[i]+U2.eta[i]-U.eta[i])/(etaM[i,kj]/100)
				# 	gradU[i,kj] <- (U2.eta[i]-U.eta[i])/(etaM[i,kj]/100)
				# }
				gradU[,kj] <- (U2.eta-U.eta)/(etaM[,kj]/100)
			}
			# 

			
				
			gradU <- gradU+(1/lambda)*(etaM - eta_prox)


			# if (u>1){
			# 	adap <- adap - gamma*(deltu + log(0.57))
			# }
			
			Z <- matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)

			for (i in 1:Dargs$NM){
				etaMc[i,] <- etaM[i,] - sigma*adap[i]*gradU[i,] + sqrt(2*sigma*adap[i])*Z[i,]
				# etaMc[i,] <- etaM[i,] + sigma*adap[i]*(gradU[i,]+1/lambda*(etaM[i,] - eta_prox[i,])) + sqrt(2*sigma*adap[i])*Z[i,]
			}
			

			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
			psiMc<-transphi(phiMc,Dargs$transform.par)
			fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
			DYF[Uargs$ind.ioM]<- -fpred
			Uc.y<-colSums(DYF)
			Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))

			#Gradient in candidate eta

			for (kj in 1:(nb.etas)){
				etaM2 <- etaMc
				phiM2 <- phiMc
				etaM2[,kj] <- etaMc[,kj] + etaMc[,kj]/100
				phiM2 <- mean.phiM[,varList$ind.eta]+etaM2
				psiM2<-transphi(phiM2,Dargs$transform.par)
				fpred2<-structural.model(psiM2, Dargs$IdM, Dargs$XM)
				DYF[Uargs$ind.ioM]<- -fpred2
				U2.y<-colSums(DYF)
				U2.eta<-0.5*rowSums(etaM2*(etaM2%*%somega))
				# for (i in 1:Dargs$NM){
				# 	gradUc[i,kj] <- (U2.eta[i]-Uc.eta[i])/(etaMc[i,kj]/100)
				# }
				
				gradUc[,kj] <- (U2.eta-Uc.eta)/(etaMc[,kj]/100)
			}

#######find the prox of g(eta_candidate)

			for(i in 1:Dargs$NM) {
			    isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			#    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-mean.phiM
			    # phii<-saemixObject["results"]["phi"][i,]
			    phii<-phiMc
			    phi1<-phii[i1.omega2]

			    phi1.opti<-optim(par=phi1, fn=proximal.function, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"],lambda=lambda,current_phi=phiMc,control = list(maxit =1))
			    # phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"],control = list(maxit = 2))
			    phi.prox[i,i1.omega2]<-phi1.opti$par
			  }

			prox.phi<-data.frame(id=id.list,phi.prox)
			prox.psi<-transphi(phi.prox,saemixObject["model"]["transform.par"])
			prox.psi<-data.frame(id=id.list,prox.psi)
			

			psic_prox <- as.matrix(prox.psi[,-c(1)])
			phic_prox <- as.matrix(prox.phi[,-c(1)])
			etac_prox <- phic_prox - mean.phiM


			
			gradUc <- gradUc + (1/lambda)*(etaMc - etac_prox)
			


			# for (i in 1:(Dargs$NM)){
			# 	propc[i,] <- ((etaMc[i,]-etaM[i,] + sigma*adap[i]*gradU[i,])/sqrt(2*sigma*adap[i]))^2
			# 	prop[i,] <- ((etaM[i,]-etaMc[i,] + sigma*adap[i]*gradUc[i,])/sqrt(2*sigma*adap[i]))^2
			# }
			
			propc <- ((etaMc-etaM + sigma*adap*gradU)/sqrt(2*sigma*adap))^2
			prop <- ((etaM-etaMc + sigma*adap*gradUc)/sqrt(2*sigma*adap))^2

			 

			P<-0.5*rowSums(prop)
			Pc<-0.5*rowSums(propc)

			
			deltu<-Uc.y-U.y+Uc.eta-U.eta + P - Pc

			ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
			
			etaM[ind,]<-etaMc[ind,]
			# for (i in 1:(nrow(phiM))) {
			# 	post_mala[[i]][u,2:(ncol(post_mala[[i]]) - 1)] <- etaM[i,]
			# }
			U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
			U.eta[ind]<-Uc.eta[ind]
			nbc2<-nbc2+length(ind)
			nt2<-nt2+Dargs$NM

			
		}
	}

	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM, post = post))
}
