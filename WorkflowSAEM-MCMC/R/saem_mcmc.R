############################### Simulation - MCMC kernels (E-step) #############################

saem.mcmc<-function(model,data,control) {
	# SAEM + MCMC
	# Input: kiter, Uargs, structural.model, mean.phi (unchanged)
	# Output: Markiv chain states
	saem <- saemix(model,data,control)
	browser()
	parpop <- tail(data.frame(saem$par), 1)
	saemixObject<-saem$object
	xinit<-saem$init

	#Options
	saemix.options<-saemixObject["options"]
	saemix.model<-saemixObject["model"]
	saemix.data<-saemixObject["data"]
	saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
	saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
	saemix.data@ntot.obs<-dim(saemix.data@data)[1]

	OLDRAND<-TRUE
  	set.seed(saemix.options$seed)
	


  	saemix.model<-xinit$saemix.model
  	Dargs<-xinit$Dargs
	saemix.model<-saemixModel(model=Dargs$structural.model,description="warfarin",type=Dargs$type
	  ,psi0=matrix(c(1,7,1,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","k"))),
	  transform.par=Dargs$transform.par,omega.init=matrix(diag(parpop[4:6]),ncol=3,byrow=TRUE),
	  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
	  byrow=TRUE))


	#Initialisation
	
	
	
	Uargs<-xinit$Uargs
	varList<-xinit$varList
	phiM<-xinit$phiM
	mean.phi<-xinit$mean.phi
	DYF<-xinit$DYF
	opt<-xinit$opt
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

	if(Dargs$type=="structural"){
		U.y<-compute.LLy_c(phiM,varList$pres,Uargs,Dargs,DYF)
	} else{
		U.y <- compute.LLy_d(phiM,Uargs,Dargs,DYF)
	}

	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	phiMc<-phiM
	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
	
	eta_list <- list(as.data.frame(matrix(nrow = saemix.options$L_mcmc,ncol = ncol(phiM))))
	for (i in 1:(nrow(phiM))) {
		eta_list[[i]] <- as.data.frame(matrix(nrow = saemix.options$L_mcmc,ncol = ncol(phiM)))
		eta_list[[i]][,ncol(eta_list[[i]])] <- i
	}

	
	if(opt$nbiter.mcmc[4]>0) {
		etaMc<-etaM
		propc <- U.eta
		prop <- U.eta
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

			fpred1<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
			gradf <- matrix(0L, nrow = length(fpred1), ncol = nb.etas) 
			for (j in 1:nb.etas) {
				psi_map2 <- psi_map
				psi_map2[,j] <- psi_map[,j]+psi_map[,j]/1000
				fpred1<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
				fpred2<-structural.model(psi_map2, Dargs$IdM, Dargs$XM)
				for (i in 1:(Dargs$NM)){
					r = which(Dargs$IdM==i)
					gradf[r,j] <- (fpred2[r] - fpred1[r])/(psi_map[i,j]/1000)
				}
			}

			gradh <- list(omega.eta,omega.eta)
			for (i in 1:Dargs$NM){
				gradh[[i]] <- gradh[[1]]
			}
			for (j in 1:nb.etas) {
				phi_map2 <- phi_map
				phi_map2[,j] <- phi_map[,j]+phi_map[,j]/1000
				psi_map2 <- transphi(phi_map2,saemixObject["model"]["transform.par"]) 
				for (i in 1:(Dargs$NM)){
					gradh[[i]][,j] <- (psi_map2[i,] - psi_map[i,])/(phi_map[i,]/1000)
				}
			}
			
			#calculation of the covariance matrix of the proposal
			Gamma <- chol.Gamma <- inv.chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
			for (i in 1:(Dargs$NM)){
				r = which(Dargs$IdM==i)
		        temp <- gradf[r,]%*%gradh[[i]]
				Gamma[[i]] <- solve(t(temp)%*%temp/(varList$pres[1])^2+solve(omega.eta))
				chol.Gamma[[i]] <- chol(Gamma[[i]])
				inv.chol.Gamma[[i]] <- solve(chol.Gamma[[i]])
				inv.Gamma[[i]] <- solve(Gamma[[i]])
			}
			
	  	} else {
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

	  	}

	  	etaM <- eta_map
	  	phiM<-etaM+mean.phiM
	  	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
	  	if(Dargs$type=="structural"){
			U.y<-compute.LLy_c(phiM,varList$pres,Uargs,Dargs,DYF)
		} else{
			U.y <- compute.LLy_d(phiM,Uargs,Dargs,DYF)
		}
	}	
	for (m in 1:saemix.options$L_mcmc) {
		if(m%%100==0){
				print(m)
		} 
	if(opt$nbiter.mcmc[4]>0) {
		for (u in 1:opt$nbiter.mcmc[4]) {
				#generate candidate eta
				for (i in 1:(Dargs$NM)){
					Mi <- rnorm(nb.etas)%*%chol.Gamma[[i]]
					etaMc[i,varList$ind.eta]<- eta_map[i,varList$ind.eta] + Mi
				}

				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc[,varList$ind.eta]
				if(Dargs$type=="structural"){
					Uc.y<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF)
				} else{
					Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
				}
				Uc.eta<-0.5*rowSums(etaMc[,varList$ind.eta]*(etaMc[,varList$ind.eta]%*%somega))

				for (i in 1:(Dargs$NM)){
					propc[i] <- 0.5*rowSums((etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%inv.Gamma[[i]])
					prop[i] <- 0.5*rowSums((etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%inv.Gamma[[i]])
				}
				
				deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind,varList$ind.eta]<-etaMc[ind,varList$ind.eta]
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]

		  	}
	}
	if(opt$nbiter.mcmc[1]>0) {
		for(u in 1:opt$nbiter.mcmc[1]) {
		# 1er noyau
			etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
			if(Dargs$type=="structural"){
				Uc.y<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF)
			} else {
				Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
			}
			deltau<-Uc.y-U.y
			ind<-which(deltau<(-1)*log(runif(Dargs$NM)))
			etaM[ind,]<-etaMc[ind,]
			U.y[ind]<-Uc.y[ind]
		}
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
					if(Dargs$type=="structural"){
						Uc.y<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF)
					} else {
						Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
					}
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
					if(Dargs$type=="structural"){
						Uc.y<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF)
					} else {
						Uc.y<-compute.LLy_d(phiMc,Uargs,Dargs,DYF)
					}
					Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
					deltu<-Uc.y-U.y+Uc.eta-U.eta
					ind<-which(deltu<(-log(runif(Dargs$NM))))
					etaM[ind,]<-etaMc[ind,]
					
			
					U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
					U.eta[ind]<-Uc.eta[ind]
					nbc2[vk2]<-nbc2[vk2]+length(ind)
					nt2[vk2]<-nt2[vk2]+Dargs$NM
				}
				
			}
			varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
		}


		for (i in 1:(nrow(phiM))) {
				eta_list[[i]][m,] <- etaM[i,]
		}
	}

	return(eta_list)
}
