################## Stochastic approximation - compute sufficient statistics (M-step) #####################
mstep_time<-function(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat) {
	# M-step - stochastic approximation
	# Input: kiter, Uargs, structural.model, DYF, phiM (unchanged)
	# Output: varList, phi, betas, suffStat (changed)
	#					mean.phi (created)
	
	# Update variances - TODO - check if here or elsewhere
	nb.etas<-length(varList$ind.eta)
	domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
	omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
	omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
	#  print(varList$omega.eta)
	chol.omega<-try(chol(omega.eta))
	d1.omega<-Uargs$LCOV[,varList$ind.eta]%*%solve(omega.eta)
	d2.omega<-d1.omega%*%t(Uargs$LCOV[,varList$ind.eta])
	comega<-Uargs$COV2*d2.omega
	
	psiM<-transphi(phiM,Dargs$transform.par)
	fpred<-structural.model(psiM, Dargs$IdM, Dargs$XM)
	
	fpred <- -fpred
	# if(Dargs$error.model=="exponential")
	# 	fpred<-log(cutoff(fpred))
	Uargs$nchains = 1
	ff<-matrix(fpred,nrow=Dargs$nobs,ncol=Uargs$nchains)
	for(k in 1:Uargs$nchains) phi[,,k]<-phiM[((k-1)*Dargs$N+1):(k*Dargs$N),]
	# overall speed similar
	#    phi<-aperm(array(phiM,c(N,nchains,3)),c(1,3,2))
	stat1<-apply(phi[,varList$ind.eta,,drop=FALSE],c(1,2),sum) # sum on columns ind.eta of phi, across 3rd dimension
	stat2<-matrix(data=0,nrow=nb.etas,ncol=nb.etas)
	stat3<-apply(phi**2,c(1,2),sum) #  sum on phi**2, across 3rd dimension
	statr<-0
	for(k in 1:Uargs$nchains) {
		phik<-phi[,varList$ind.eta,k]
		stat2<-stat2+t(phik)%*%phik
		fk<-ff[,k]
	}
	# Update sufficient statistics
	suffStat$statphi1<-suffStat$statphi1+opt$stepsize[kiter]*(stat1/Uargs$nchains-suffStat$statphi1)
	suffStat$statphi2<-suffStat$statphi2+opt$stepsize[kiter]*(stat2/Uargs$nchains-suffStat$statphi2)
	suffStat$statphi3<-suffStat$statphi3+opt$stepsize[kiter]*(stat3/Uargs$nchains-suffStat$statphi3)
	suffStat$statrese<-suffStat$statrese+opt$stepsize[kiter]*(statr/Uargs$nchains-suffStat$statrese)
	
	############# Maximisation
	##### fixed effects
	if (opt$flag.fmin && kiter>=0) {
		temp<-d1.omega[Uargs$ind.fix11,]*(t(Uargs$COV1)%*%(suffStat$statphi1-Uargs$dstatCOV[,varList$ind.eta]))
		betas[Uargs$ind.fix11]<-solve(comega[Uargs$ind.fix11,Uargs$ind.fix11],rowSums(temp)) 
		# ECO TODO: utiliser optimise dans le cas de la dimension 1
		
		beta0<-optim(par=betas[Uargs$ind.fix10],fn=compute.Uy_time,phiM=phiM,args=Uargs,Dargs=Dargs,DYF=DYF,control=list(maxit=opt$maxim.maxiter))$par # else
		betas[Uargs$ind.fix10]<-betas[Uargs$ind.fix10]+opt$stepsize[kiter]*(beta0-betas[Uargs$ind.fix10])
	} else {
		temp<-d1.omega[Uargs$ind.fix1,]*(t(Uargs$COV1)%*%(suffStat$statphi1-Uargs$dstatCOV[,varList$ind.eta]))
		betas[Uargs$ind.fix1]<-solve(comega[Uargs$ind.fix1,Uargs$ind.fix1],rowSums(temp)) 
	}

	varList$MCOV[Uargs$j.covariate]<-betas
	mean.phi<-Uargs$COV %*% varList$MCOV
	e1.phi<-mean.phi[,varList$ind.eta,drop=FALSE]
	
	# Covariance of the random effects
	omega.full<-matrix(data=0,nrow=Uargs$nb.parameters,ncol=Uargs$nb.parameters)
	omega.full[varList$ind.eta,varList$ind.eta]<-suffStat$statphi2/Dargs$N + t(e1.phi)%*%e1.phi/Dargs$N - t(suffStat$statphi1)%*%e1.phi/Dargs$N - t(e1.phi)%*%suffStat$statphi1/Dargs$N
	varList$omega[Uargs$indest.omega]<-omega.full[Uargs$indest.omega]
	

	# Simulated annealing (applied to the diagonal elements of omega)
	if (kiter<=0) {
		diag.omega.full<-mydiag(omega.full)
		vec1<-diag.omega.full[Uargs$i1.omega2]
		vec2<-varList$diag.omega[Uargs$i1.omega2]*opt$alpha1.sa
		idx<-as.integer(vec1<vec2)
		varList$diag.omega[Uargs$i1.omega2]<-idx*vec2+(1-idx)*vec1
		varList$diag.omega[Uargs$i0.omega2]<-varList$diag.omega[Uargs$i0.omega2]*opt$alpha0.sa
	} else {
		varList$diag.omega<-mydiag(varList$omega)
	}
	varList$omega<-varList$omega-mydiag(mydiag(varList$omega))+mydiag(varList$diag.omega)

	return(list(varList=varList,mean.phi=mean.phi,phi=phi,betas=betas,suffStat=suffStat))
}
