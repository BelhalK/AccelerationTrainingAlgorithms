############################### Simulation - MCMC kernels (E-step) #############################

estep_incremental<-function(kiter, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM,saemixObject,l,ind_rand) {
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
	saemix.options<-saemixObject["options"]
	map_range <- saemix.options$map.range


	saemix.options<-saemixObject["options"]
  	saemix.model<-saemixObject["model"]
  	saemix.data<-saemixObject["data"]
  	saemix.options$map <- TRUE
  	saemixObject["results"]["omega"] <- omega.eta
  	saemixObject["results"]["mean.phi"] <- mean.phi
  	saemixObject["results"]["phi"] <- phiM
  	i1.omega2<-varList$ind.eta
    iomega.phi1<-solve(saemixObject["results"]["omega"][i1.omega2,i1.omega2])
  	id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  	xind<-saemixObject["data"]["data"][,saemixObject["data"]["name.predictors"], drop=FALSE]
  	yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  	id.list<-unique(id)
  	phi.map<-saemixObject["results"]["mean.phi"]
  	eta_map <- phi.map
  	indchosen <- l[ind_rand]
# Sampling strategy (MAP calculation)
if (kiter <= 0){ #if rwm
  	# if (kiter <= length(map_range) && length(ind_rand)!=Dargs$NM){
	 for(i in 1:saemixObject["data"]["N"]) {
	    isuj<-id.list[i]
	    xi<-xind[id==isuj,,drop=FALSE]
	    yi<-yobs[id==isuj]
	    idi<-rep(1,length(yi))
	    mean.phi1<-mean.phiM[i,i1.omega2]
	    phii<-saemixObject["results"]["phi"][i,]
	    phi1<-phii[i1.omega2]
	    phi1.opti<-optim(par=phi1, fn=conditional.distribution_c, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=varList$pres, err=saemixObject["model"]["error.model"])
	    phi.map[i,i1.omega2]<-phi1.opti$par
	}
	#rep the map nchains time
	phi.map <- phi.map[rep(seq_len(nrow(phi.map)),Uargs$nchains ), ]

  	map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
	map.psi<-data.frame(id=id.list,map.psi)
	map.phi<-data.frame(id=id.list,phi.map)
	psi_map <- as.matrix(map.psi[,-c(1)])
	phi_map <- as.matrix(map.phi[,-c(1)])
	eta_map <- phi_map - mean.phiM

	weight <- eta_map[,1]
	gamma = saemix.options$gamma
	
	for (m in 1:Dargs$NM){
		weight[m] <- exp(gamma*eta_map[m,2]^2)
		# weight[m] <- exp(gamma*eta_map[m,2])
	}
	weight <- weight/sum(weight)
	nb.replacement <- length(ind_rand)
	indchosen <- sample(1:Dargs$NM, size = nb.replacement, replace = FALSE, prob = weight)
	# indchosen <- sample(1:Dargs$NM, size = nb.replacement, replace = FALSE)
	block <- setdiff(1:Dargs$NM, indchosen)
	
	# print(kiter)
	# print(indchosen)
	
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	etaM[indchosen,] <- eta_map[indchosen,] #if rwm
	phiM <- etaM + mean.phiM
	phiMc<-phiM

} else {
	
	#indchosen <- 1:Dargs$NM
	block <- NULL
	for (m in 1:Uargs$nchains){	
		block <- list.append(block,setdiff(1:Dargs$N, l[ind_rand])+(m-1)*Dargs$N)
	}
	chosen <- NULL
	for (m in 1:Uargs$nchains){	
		chosen <- list.append(chosen, l[ind_rand]+(m-1)*Dargs$N)
	}
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	phiMc<-phiM

	chosen <- unique(chosen)
}

if(Dargs$type=="structural"){
		U.y<-compute.LLy_c(phiM,varList$pres,Uargs,Dargs,DYF,1:Dargs$N)
	} else{
		U.y <- compute.LLy_d(phiM,Uargs,Dargs,DYF)
}

	# block <- setdiff(1:Dargs$NM, l[ind_rand])	
	
	if (!(kiter %in% map_range)){
		for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau
			etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
			if(Dargs$type=="structural"){
				Uc.y<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF,chosen)
			} else {
				Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
			}
			deltau<-Uc.y-U.y
			deltau[block] = 1000000
			ind<-which(deltau<(-1)*log(runif(Dargs$NM)))
			# print(length(ind)/length(indchosen))
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
					# psiMc<-transphi(phiMc,Dargs$transform.par)
					if(Dargs$type=="structural"){
						Uc.y<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF,chosen)
					} else {
						Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
					}
					Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
					deltu<-Uc.y-U.y+Uc.eta-U.eta
					deltu[block] = 1000000
					ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
					# print(length(ind)/length(indchosen))
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
					# psiMc<-transphi(phiMc,Dargs$transform.par)
					if(Dargs$type=="structural"){
						Uc.y<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF,chosen)
					} else {
						Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
					}
					Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
					deltu<-Uc.y-U.y+Uc.eta-U.eta
					deltu[block] = 1000000
					ind<-which(deltu<(-log(runif(Dargs$NM))))
					# print(length(ind)/length(indchosen))
					etaM[ind,]<-etaMc[ind,]
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

	if(opt$nbiter.mcmc[4]>0 & kiter %in% map_range){
		saemix.options<-saemixObject["options"]
	  	saemix.model<-saemixObject["model"]
	  	saemix.data<-saemixObject["data"]
	  	saemix.options$map <- TRUE
	  	saemixObject["results"]["omega"] <- omega.eta
	  	saemixObject["results"]["mean.phi"] <- mean.phi
	  	saemixObject["results"]["phi"] <- phiM
	  	i1.omega2<-varList$ind.eta
	    iomega.phi1<-solve(saemixObject["results"]["omega"][i1.omega2,i1.omega2])
	  	id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
	  	xind<-saemixObject["data"]["data"][,saemixObject["data"]["name.predictors"], drop=FALSE]
	  	yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
	  	id.list<-unique(id)
	  	phi.map<-saemixObject["results"]["mean.phi"]

		if(Dargs$type=="structural"){
			#MAP calculation
		 	# for(i in 1:saemixObject["data"]["N"]) {
			# phiM <- saemixObject["results"]["phi"]
			# psiM <- transphi(phiM,saemixObject["model"]["transform.par"])
			
			# fpred1 <- structural.model(psiM,Dargs$IdM,Dargs$XM)
			# tempsiM <- cbind(unique(Dargs$IdM), psiM)
			# colnames(tempsiM) <- c("id",colnames(omega.eta))
			# fpred2 <- computePredictions(data.frame(tempsiM))[[1]]
			
			# structural.model(psiM[10,],rep(1,11),Dargs$XM[Dargs$IdM==10,])
			# computePredictions(data.frame(tempsiM)[10,], individualIds=10)[[1]]
			# computePredictions(data.frame(tempsiM)[26,], individualIds=26)[[1]]
			# computePredictions(data.frame(tempsiM)[28,], individualIds=28)[[1]]
			# structural.model(psiM,Dargs$IdM[Dargs$IdM==28],Dargs$XM[which(Dargs$IdM==28),])
			for(i in chosen) {
			    isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-mean.phiM[i,i1.omega2]
			    phii<-saemixObject["results"]["phi"][i,]
			    phi1<-phii[i1.omega2]
			    phi1.opti<-optim(par=phi1, fn=conditional.distribution_c, 
			    	phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,
			    	iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], 
			    	model=saemixObject["model"]["model"], pres=varList$pres, 
			    	err=saemixObject["model"]["error.model"], index.indiv=i)
			    phi.map[i,i1.omega2]<-phi1.opti$par
			}
			phimap1<- phi.map
			# for(i in chosen) {
			#     isuj<-id.list[i]
			#     xi<-xind[id==isuj,,drop=FALSE]
			#     yi<-yobs[id==isuj]
			#     idi<-rep(1,length(yi))
			#     mean.phi1<-mean.phiM[i,i1.omega2]
			#     phii<-saemixObject["results"]["phi"][i,]
			#     phi1<-phii[i1.omega2]
			#     phi1.opti<-optim(par=phi1, fn=conditional.distribution_c_test, 
			#     	phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,
			#     	iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], 
			#     	model=saemixObject["model"]["model"], pres=varList$pres, 
			#     	err=saemixObject["model"]["error.model"])
			#     phi.map[i,i1.omega2]<-phi1.opti$par
			# }
			# phimap2<- phi.map
			# phimap2 - phi.map
			#rep the map nchains time
			phi.map <- phi.map[rep(seq_len(nrow(phi.map)),Uargs$nchains ), ]
		  	map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
			map.psi<-data.frame(id=id.list,map.psi)
			map.phi<-data.frame(id=id.list,phi.map)
			psi_map <- as.matrix(map.psi[,-c(1)])
			phi_map <- as.matrix(map.phi[,-c(1)])
			eta_map <- phi_map - mean.phiM
		
			#gradient at the map estimation
			tempsi_map <- cbind(unique(Dargs$IdM), psi_map)
			colnames(tempsi_map) <- c("id",colnames(omega.eta))
			fpred1<-computePredictions(data.frame(tempsi_map))[[1]]

			# fpred1<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
			gradf <- matrix(0L, nrow = length(fpred1), ncol = nb.etas) 
			r <- which(Dargs$IdM %in% chosen)
			for (j in 1:nb.etas) {
				phi_map2 <- phi_map
				phi_map2[,j] <- phi_map[,j]+1/100;
				psi_map2 <- transphi(phi_map2,saemixObject["model"]["transform.par"]) 
				tempsi_map <- cbind(unique(Dargs$IdM), psi_map)
				tempsi_map2 <- cbind(unique(Dargs$IdM), psi_map2)
				colnames(tempsi_map) <- c("id",colnames(omega.eta))
				colnames(tempsi_map2) <- c("id",colnames(omega.eta))
				fpred1 <- fpred2 <- Dargs$yM
				# if (kiter==2) browser()
				fpred1[r]<-computePredictions(data.frame(tempsi_map)[chosen,], individualIds=chosen)[[1]]
				fpred2[r]<-computePredictions(data.frame(tempsi_map2)[chosen,], individualIds=chosen)[[1]]

				for (i in chosen){
					gradf[which(Dargs$IdM == i),j] <- (fpred2[which(Dargs$IdM == i)] - fpred1[which(Dargs$IdM == i)])/(1/10)
				}
			}
			#calculation of the covariance matrix of the proposal
			Gamma <- list(omega.eta,omega.eta)
			z <- matrix(0L, nrow = length(fpred1), ncol = 1) 
			for (i in 1:Dargs$NM){
				Gamma[[i]] <- omega.eta
			}

			for (i in chosen){
				Gamma[[i]] <- solve(t(gradf[which(Dargs$IdM == i),])%*%gradf[which(Dargs$IdM == i),]/(varList$pres[1])^2+solve(omega.eta))
				# Gamma[[i]] <- omega.eta
			}

			etaM[chosen,] <- eta_map[chosen,]
			for (u in 1:opt$nbiter.mcmc[4]) {
				etaMc<-etaM
				propc <- U.eta
				prop <- U.eta
				#generate candidate eta
				for (i in chosen){
					M <- matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)%*%chol(Gamma[[i]])
					etaMc[i,varList$ind.eta]<- eta_map[i,varList$ind.eta] +M[i,]
				}

				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc[,varList$ind.eta]
				Uc.y<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF,chosen)
				Uc.eta<-0.5*rowSums(etaMc[,varList$ind.eta]*(etaMc[,varList$ind.eta]%*%somega))

				# for (i in 1:(Dargs$NM)){
				for (i in chosen){
					propc[i] <- 0.5*rowSums((etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%solve(Gamma[[i]]))
					prop[i] <- 0.5*rowSums((etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%solve(Gamma[[i]]))
				}
				deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
				deltu[block] = 1000000
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind,varList$ind.eta]<-etaMc[ind,varList$ind.eta]
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
			}
		} else{
			# MAP calculation
		  	for(i in 1:saemixObject["data"]["N"]) {
			    isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			#    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-mean.phiM[i,i1.omega2]
			    phii<-saemixObject["results"]["phi"][i,]
			    phi1<-phii[i1.omega2]
			    phi1.opti<-optim(par=phi1, fn=conditional.distribution_d, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"])
			    # phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"],control = list(maxit = 2))
			    phi.map[i,i1.omega2]<-phi1.opti$par
			}
			#rep the map nchains time
			phi.map <- phi.map[rep(seq_len(nrow(phi.map)),Uargs$nchains ), ] 
		  	map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
			map.psi<-data.frame(id=id.list,map.psi)
			map.phi<-data.frame(id=id.list,phi.map)

			psi_map <- as.matrix(map.psi[,-c(1)])
			phi_map <- as.matrix(map.phi[,-c(1)])
			eta_map <- phi_map[,varList$ind.eta] - mean.phiM[,varList$ind.eta]
			
			#gradient at the map estimation
			gradp <- matrix(0L, nrow = Dargs$NM, ncol = nb.etas) 

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
					gradp[i,j] <- (l2[i] - l1[i])/(phi_map[i,j]/100)
				}
			}

			#calculation of the covariance matrix of the proposal
			fpred<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
			DYF[Uargs$ind.ioM]<- fpred
			denom <- colSums(DYF)
			
			Gamma <- chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
			z <- matrix(0L, nrow = length(fpred), ncol = 1) 
			for (i in 1:(Dargs$NM)){
				Gamma[[i]] <- solve(gradp[i,]%*%t(gradp[i,])/denom[i]^2+solve(omega.eta))
				chol.Gamma[[i]] <- chol(Gamma[[i]])
				inv.Gamma[[i]] <- solve(Gamma[[i]])
			}
			
			etaM[chosen,] <- eta_map[chosen,]
			for (u in 1:opt$nbiter.mcmc[4]) {
				etaMc<-etaM
				propc <- U.eta
				prop <- U.eta
				#generate candidate eta
				for (i in 1:(Dargs$NM)){
					Mi <- rnorm(nb.etas)%*%chol.Gamma[[i]]
					etaMc[i,]<- eta_map[i,] +Mi
				}
				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc

				Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				
				for (i in 1:(Dargs$NM)){
					propc[i] <- 0.5*rowSums((etaMc[i,]-eta_map[i,])*(etaMc[i,]-eta_map[i,])%*%inv.Gamma[[i]])
					prop[i] <- 0.5*rowSums((etaM[i,]-eta_map[i,])*(etaM[i,]-eta_map[i,])%*%inv.Gamma[[i]])
				}

				deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
				deltu[block] = 1000000
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind]<-etaMc[ind]
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
			}
		}
	}
	
	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM[,varList$ind.eta]
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM, map = eta_map, indchosen = indchosen))
}
