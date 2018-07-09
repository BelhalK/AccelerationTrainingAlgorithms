############################### Simulation - MCMC kernels (E-step) #############################

indiv.variational.inference<-function(model,data,control=list()) {
	# E-step - simulate unknown parameters
	# Input: kiter, Uargs, structural.model, mean.phi (unchanged)
	# Output: varList, DYF, phiM (changed)
	kiter <- 1
	saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
	saemix.options<-saemixObject["options"]
  	saemix.model<-saemixObject["model"]
  	saemix.data<-saemixObject["data"]
	# Initialising random generator
	OLDRAND<-TRUE
	set.seed(saemix.options$seed)
	#intitialisation
	# xinit<-initialiseMainAlgo(saemix.data,saemix.model,saemix.options)
	xinit<-initialiseMainAlgo(saemix.data,saemix.model,saemix.options)
	saemix.model<-xinit$saemix.model
	Dargs<-xinit$Dargs
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
	somega <- solve(omega.eta)
	# "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)

	VK<-rep(c(1:nb.etas),2)
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]
	saemix.options<-saemixObject["options"]

	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	phiMc<-phiM

	etaMc<-etaM
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
	Gammamap <- list(omega.eta,omega.eta)
	for (i in 1:(Dargs$NM)){
		r = which(Dargs$IdM==i)
        temp <- gradf[r,]%*%gradh[[i]]
		Gammamap[[i]] <- solve(t(temp)%*%temp/(varList$pres[1])^2+solve(omega.eta))
	}

##WARFA
	indiv <- control$indiv.index
	obs <- Dargs$yM[Dargs$IdM==indiv]
	dose <- unique(Dargs$XM[Dargs$IdM==indiv,1])
	time <- Dargs$XM[Dargs$IdM==indiv,2]
	mean.psiM <- transphi(mean.phiM,Dargs$transform.par)
	stan.model <- control$modelstan
	stan_data <- list(N = length(obs),concentration = obs
					,time = time, dose = dose,
					beta1_pop=mean.phiM[indiv,1],beta2_pop=mean.phiM[indiv,2],beta3_pop=mean.phiM[indiv,3],
					omega_beta1=sqrt(omega.eta[1,1]),omega_beta2=sqrt(omega.eta[2,2]),omega_beta3=sqrt(omega.eta[3,3]),
					pres=sqrt(varList$pres[1]))
	fit <- vb(stan.model, data = stan_data, iter = 100000)
	fit_samples = extract(fit)

	psiMstan <- tail(fit_samples$beta,L_mcmc)
	phiMstan<-transpsi(psiMstan,Dargs$transform.par)
	etaMstan <- phiMstan
	etaMstan[,1] <- phiMstan[,1] - mean.phiM[i,1]
	etaMstan[,2] <- phiMstan[,2] - mean.phiM[i,2]
	etaMstan[,3] <- phiMstan[,3] - mean.phiM[i,3]
	eta_map[indiv,]
	colMeans(etaMstan)
	mu <- colMeans(etaMstan)
	Gamma <- cov(etaMstan)

}

###LINEAR
	# indiv <- control$indiv.index
	# obs <- Dargs$yM[Dargs$IdM==indiv]
	# design <- as.data.frame(matrix(0, ncol = ncol(etaM), nrow = length(obs)))
	# design[,1] <- 1
	# design[,2] <- Dargs$XM[Dargs$IdM==indiv,]
	# design <- as.matrix(design)
	# mean.psiM <- transphi(mean.phiM,Dargs$transform.par)
	# stan.model <- control$modelstan
	# stan_data <- list(N = length(obs),height = obs
	# 				,age = design[,2],
	# 				beta1_pop=mean.phiM[indiv,1],beta2_pop=mean.phiM[indiv,2],
	# 				omega_beta1=omega.eta[1,1],omega_beta2=omega.eta[2,2],
	# 				pres=sqrt(varList$pres[1]))
	# # fit <- sampling(stan.model, data = stan_data)
	# fit <- vb(stan.model, data = stan_data, iter = 50000)
	# fit_samples = extract(fit)
	# psiMstan <- tail(fit_samples$beta,L_mcmc)
	# phiMstan<-transpsi(psiMstan,Dargs$transform.par)
	# etaMstan <- phiMstan
	# etaMstan[,1] <- phiMstan[,1] - mean.phiM[,1]
	# etaMstan[,2] <- phiMstan[,2] - mean.phiM[,2]



#### RTTE
if(Dargs$type=="likelihood"){
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
	
	Gammamap <- list(omega.eta,omega.eta)
	z <- matrix(0L, nrow = length(fpred), ncol = 1) 
	for (i in 1:(Dargs$NM)){
		Gammamap[[i]] <- solve(gradp[i,]%*%t(gradp[i,])/denom[i]^2+solve(omega.eta))
	}

	indiv <- control$indiv.index
	stan.model <- control$modelstan
	T <- Dargs$XM[Dargs$IdM==indiv,1]
	T_c <- 20
	event_times <- T[!(T %in% c(0, T_c))]
	cens_times <- T[T == T_c]
	N_e <- length(event_times)
	N_c <- length(cens_times)
	mean.psiM <- transphi(mean.phiM,Dargs$transform.par)
	stan_data <- list(N_e = N_e, N_c = N_c
					,event_times = event_times, cens_times = cens_times,
					alpha_pop=mean.phiM[indiv,1],sigma_pop=mean.phiM[indiv,2],
					omega_alpha=sqrt(omega.eta[1,1]),omega_sigma=sqrt(omega.eta[2,2]))

	fit <- vb(stan.model, data = stan_data)
	fit_samples = extract(fit)
	psiMstan <- tail(fit_samples$beta,L_mcmc)
	phiMstan<-transpsi(psiMstan,Dargs$transform.par)
	etaMstan <- phiMstan
	etaMstan[,1] <- phiMstan[,1] - mean.phiM[indiv,1]
	etaMstan[,2] <- phiMstan[,2] - mean.phiM[indiv,2]
	print(eta_map[indiv,])
	print(colMeans(etaMstan))
	colMeans(phiMstan)
	browser()
	mu <- colMeans(etaMstan)
	Gamma <- cov(etaMstan)
}

	return(list(mu=mu, Gamma = Gamma, map = eta_map, Gammamap = Gammamap))
}
