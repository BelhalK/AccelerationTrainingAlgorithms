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


	#Initialization
	L <- 100 #nb iterations MONTE CARLO
	rho <- 10e-10 #gradient ascent stepsize
	K <- control$nbiter.gd
	i <- 10
	#VI to find the right mean mu (gradient descent along the elbo)
	# mu <- list(etaM[i,],etaM[i,])
	
	# mu[[1]][1] = 169
	# mu[[1]][2] = 1
	# #if Gamma fixed to Laplace Gamma
	# trueGamma <- control$Gamma.laplace
	# chol.Gamma <- chol(trueGamma[[i]])
	# inv.Gamma <- solve(trueGamma[[i]])

	# obs <- Dargs$yM[Dargs$IdM==i]
	# design <- as.data.frame(matrix(0, ncol = ncol(etaM), nrow = length(obs)))
	# design[,1] <- 1
	# design[,2] <- Dargs$XM[Dargs$IdM==i,]
	# design <- as.matrix(design)

	# for (k in 1:K) {
	# 	if (k%%10==0) print(k)
	# 		sample <- list(etaM[i,],etaM[i,])  #list of samples for monte carlo integration
	# 		estim <- list(etaM[i,],etaM[i,]) #noisy gradient of the ELBO
	# 		gradlogq <- etaM[i,]
	# 		for (l in 1:L) {
	# 			sample[[l]] <- mu[[k]] +matrix(rnorm(nb.etas), ncol=nb.etas)%*%chol.Gamma
	# 			phiMc[i,varList$ind.eta]<-mean.phiM[i,varList$ind.eta]+sample[[l]]
	# 			psiMc<-transphi(phiMc,Dargs$transform.par)
	# 			fpred<-structural.model(psiMc, Dargs$IdM, Dargs$XM)
	# 			if(Dargs$error.model=="exponential")
	# 				fpred<-log(cutoff(fpred))
	# 			gpred<-error(fpred,varList$pres)
	# 			DYF[Uargs$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
	# 			#Log complete computation
	# 			sumDYF <- colSums(DYF)
	# 			logp <- sumDYF[i] + 0.5*rowSums(sample[[l]]*(sample[[l]]%*%somega))
	# 			#Log proposal computation
	# 			logq <- 0.5*rowSums((sample[[l]] - mu[[k]])*((sample[[l]] - mu[[k]])%*%inv.Gamma))
	# 			#gradlogq computation
	# 			for (j in 1:nb.etas) {
	# 				mu2 <- mu[[k]]
	# 				mu2[j] <- mu[[k]][j] + mu[[k]][j]/100
	# 				logq2 <- 0.5*rowSums((sample[[l]] - mu2)*((sample[[l]] - mu2)%*%inv.Gamma))
	# 				gradlogq[j] <- (logq2 - logq) / (mu[[k]][j]/100)
	# 			}
	# 			estim[[l]] <- (logq - logp)*gradlogq
	# 		}
	# 		grad_mu_elbo <- 1/L*Reduce("+", estim) 
			
	# 		#Gradient ascent along that gradient
	# 		mu[[k+1]] <- mu[[k]] - rho*grad_mu_elbo
	# }

	# k <- 1
	# first <- t(omega.eta%*%(mu[[k]]- mean.phiM[i,]))
	# second <- (-t(obs)%*%design - mu[[k]]%*%t(design)%*%design)/varList$pres[1]
	# grad_mu_elbo <- first + second
			
	# #Gradient ascent along that gradient
	# mu[[k+1]] <- mu[[k]] - rho*grad_mu_elbo
	# for (k in 2:K) {
	# 	if (k%%10==0) print(k)
	# 		first <- t(omega.eta%*%(t(mu[[k]])- mean.phiM[i,]))
	# 		second <- (-t(obs)%*%design - mu[[k]]%*%t(design)%*%design)/varList$ pres[1]
	# 		grad_mu_elbo <- first + second
			
	# 		#Gradient ascent along that gradient
	# 		mu[[k+1]] <- mu[[k]] - rho*grad_mu_elbo
	# }

###LINEAR
	# i <- 10
	# obs <- Dargs$yM[Dargs$IdM==i]
	# design <- as.data.frame(matrix(0, ncol = ncol(etaM), nrow = length(obs)))
	# design[,1] <- 1
	# design[,2] <- Dargs$XM[Dargs$IdM==i,]
	# design <- as.matrix(design)
	
	# stan.model <- control$modelstan
	# stan_data <- list(N = length(obs),height = obs
	# 				,age = design[,2],
	# 				beta1_pop=mean.phiM[i,1],beta2_pop=mean.phiM[i,2],
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

###WARFA
	i <- 10
	obs <- Dargs$yM[Dargs$IdM==i]
	design <- as.data.frame(matrix(0, ncol = ncol(etaM), nrow = length(obs)))
	design[,1] <- Dargs$XM[1,1]
	design[,2] <- Dargs$XM[Dargs$IdM==i,2]
	design <- as.matrix(design)
	
	mean.psiM <- transphi(mean.phiM,saemixObject["model"]["transform.par"])
	stan.model <- control$modelstan
	obs <- obs[-1]
	stan_data <- list(N = length(obs),concentration = obs
					,time = design[-1,2],
					beta1_pop=mean.psiM[i,1],beta2_pop=mean.psiM[i,2],beta3_pop=mean.psiM[i,3],
					omega_beta1=omega.eta[1,1],omega_beta2=omega.eta[2,2],omega_beta3=omega.eta[3,3],
					pres=sqrt(varList$pres[1]))


	fit <- vb(stan.model, data = stan_data, iter = 100000)
	fit_samples = extract(fit)
	
	psiMstan <- tail(fit_samples$beta,L_mcmc)
	phiMstan<-transpsi(psiMstan,Dargs$transform.par)
	etaMstan <- phiMstan
	etaMstan[,1] <- phiMstan[,1] - mean.phiM[i,1]
	etaMstan[,2] <- phiMstan[,2] - mean.phiM[i,2]
	etaMstan[,3] <- phiMstan[,3] - mean.phiM[i,3]
	
	eta_map[i,]
	colMeans(etaMstan)
	mu <- colMeans(etaMstan)
	Gamma <- cov(etaMstan)
	
	return(list(mu=mu, Gamma = Gamma))
}
