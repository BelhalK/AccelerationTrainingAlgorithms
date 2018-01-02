############################### Simulation - MCMC kernels (E-step) #############################

estep<-function(kiter, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM,saemixObject) {
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

	if(opt$flag.fmin){
		omega0.eta <- varList$omega[varList$ind0.eta,varList$ind0.eta,drop=FALSE]
		chol.omega0<-try(chol(omega0.eta))
		somega0<-solve(omega0.eta)
	}
	
	# "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)
	
	VK<-rep(c(1:nb.etas),2)
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
	saemix.options<-saemixObject["options"]
	map_range <- saemix.options$map.range

	if(Dargs$type=="structural"){
		U.y<-compute.LLy_c(phiM,varList$pres,Uargs,Dargs,DYF)
	} else{
		U.y <- compute.LLy_d(phiM,Uargs,Dargs,DYF)
	}
	
	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	etaM0<-phiM[1,varList$ind0.eta]-mean.phiM[1,varList$ind0.eta,drop=FALSE]
	phiMc<-phiM
	if (!(kiter %in% map_range)){
	for(u in 1:opt$nbiter.mcmc[1]) { # 1er noyau
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
#MCMC on the fixed param (only 1 and not N)
U0.y <- sum(U.y)
		

if(kiter < opt$nbiter.saemix[1]/2){
if(opt$flag.fmin & opt$nbiter.mcmc[4]>0){
	for(u in 1:opt$nbiter.mcmc[4]) { 
		# browser()
			etaMc0 <- matrix(rnorm(length(varList$ind0.eta)),ncol=length(varList$ind0.eta))%*%chol.omega0
			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
			phiMc[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]+etaMc0
			if(Dargs$type=="structural"){
				Uc0.y<-sum(compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF))
			} else {
				Uc0.y<-sum(compute.LLy_d(phiMc,Uargs,Dargs,DYF))
			}
			delt<-Uc0.y-U0.y
			ind<-which(delt<(-1)*log(runif(length(varList$ind0.eta))))
			print(ind)
			etaM0[ind,]<-etaMc0[ind,]
			U0.y[ind]<-Uc0.y[ind]
		}

	# U0.eta<-0.5*rowSums(etaM0*(etaM0%*%somega0))
	# nt2<-nbc2<-matrix(data=0,nrow=length(varList$ind0.eta),ncol=1)
	# nrs2<-1
	# for (u in 1:opt$nbiter.mcmc[4]) {
	# 	for(vk2 in 1:length(varList$ind0.eta)) {
	# 		etaMc0<-etaM0
	# 		#				cat('vk2=',vk2,' nrs2=',nrs2,"\n")
	# 		etaMc0[,vk2]<-etaM0[,vk2]+matrix(rnorm(length(varList$ind0.eta)), ncol=nrs2)%*%chol.omega0 
	# 		phiMc[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]+etaMc0
	# 		psiMc<-transphi(phiMc,Dargs$transform.par)
	# 		if(Dargs$type=="structural"){
	# 			Uc0.y<-sum(compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF))
	# 		} else {
	# 			Uc0.y<-sum(compute.LLy_d(phiMc,Uargs,Dargs,DYF))
	# 		}
	# 		Uc0.eta<-0.5*rowSums(etaMc0*(etaMc0%*%somega0))
	# 		deltu<-Uc0.y-U0.y+Uc0.eta-U0.eta
	# 		ind<-which(deltu<(-1)*log(runif(length(varList$ind0.eta))))
	# 		print(ind)
	# 		etaM0[ind,]<-etaMc0[ind,]
	# 		U0.y[ind]<-Uc0.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
	# 		U0.eta[ind]<-Uc0.eta[ind]
	# 		nbc2[vk2]<-nbc2[vk2]+length(ind)
	# 	}
	# }
}
}



if(Dargs$type=="structural"){
	if(opt$nbiter.mcmc[4]>0 & kiter %in% map_range){
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

		  	# i1.omega2<-saemixObject["model"]["indx.omega"]
		  	i1.omega2<-varList$ind.eta
		    iomega.phi1<-solve(saemixObject["results"]["omega"][i1.omega2,i1.omega2])
		  	id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
		  	xind<-saemixObject["data"]["data"][,saemixObject["data"]["name.predictors"], drop=FALSE]
		  	yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
		  	id.list<-unique(id)
		  	phi.map<-saemixObject["results"]["phi"]


			K1 <- saemix.options$nbiter.saemix[1]
			K2 <- saemix.options$nbiter.saemix[2]
			K_gd <- 1
			gd_step <- saemix.options$step.gd
			if (kiter %in% map_range){
			  	for(i in 1:Dargs$NM) {
				    isuj<-id.list[i]
				    xi<-xind[id==isuj,,drop=FALSE]
				    yi<-yobs[id==isuj]
				    idi<-rep(1,length(yi))
				    mean.phi1<-saemixObject["results"]["mean.phi"][i,i1.omega2]
				    phii<-saemixObject["results"]["phi"][i,]
				    phi1<-phii[i1.omega2]
				    phi1.opti<-optim(par=phi1, fn=conditional.distribution_c, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"])
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
			fpred1<-structural.model(psi_map, Dargs$IdM, Dargs$XM)
			gradf <- matrix(0L, nrow = length(fpred1), ncol = nb.etas) 


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
			z <- matrix(0L, nrow = length(fpred1), ncol = 1) 
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
				etaMc<-etaM
				etaM <- eta_map
				propc <- U.eta
				prop <- U.eta
				#generate candidate eta
					
				for (i in 1:(Dargs$NM)){
					M <- matrix(rnorm(Dargs$NM*nb.etas), ncol=nb.etas)%*%chol(Gamma[[i]])
					etaMc[i,varList$ind.eta]<- eta_map[i,varList$ind.eta] +M[i,]
				}

				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc[,varList$ind.eta]
				Uc.y<-compute.LLy_c(phiMc,varList$pres,Uargs,Dargs,DYF)
				Uc.eta<-0.5*rowSums(etaMc[,varList$ind.eta]*(etaMc[,varList$ind.eta]%*%somega))


				for (i in 1:(Dargs$NM)){
					propc[i] <- 0.5*rowSums((etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%solve(Gamma[[i]]))
					prop[i] <- 0.5*rowSums((etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%solve(Gamma[[i]]))
				}


				deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))

				etaM[ind,varList$ind.eta]<-etaMc[ind,varList$ind.eta]
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
			}
		}
} else{
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
		i1.omega2<-varList$ind.eta
	  	# i1.omega2<-saemixObject["model"]["indx.omega"]
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
		    phi1.opti<-optim(par=phi1, fn=conditional.distribution_d, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"])
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
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				for (i in 1:(Dargs$NM)){
					propc[i] <- 0.5*rowSums((etaMc[i,]-eta_map[i,])*(etaMc[i,]-eta_map[i,])%*%inv.Gamma[[i]])
					prop[i] <- 0.5*rowSums((etaM[i,]-eta_map[i,])*(etaM[i,]-eta_map[i,])%*%inv.Gamma[[i]])
				}

				deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind]<-etaMc[ind]
				U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM


			}
		}
	}

}
	
	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM[,varList$ind.eta]
	if(opt$flag.fmin){
		phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]+etaM0
	}
	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM))
}
