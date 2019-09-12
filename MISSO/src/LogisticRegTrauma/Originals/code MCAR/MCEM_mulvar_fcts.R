#                                   Functions


#----------------------MCEM with Acceptation-rejection sampling-----------------------------
#mcem_ar=function(X,Y,X.cc,beta.start,mu.start,sig2.start,maxruns,tol_em,print_iter=TRUE){
mcem_ar=function(X,Y, maxruns=1000,tol_em=1e-5,print_iter=TRUE, NBITERINIT= 20, mc.sizefin=1000){
  p=dim(X)[2]
  if(any(apply(is.na(X),1,sum)==p)){
    i_allNA=which(apply(is.na(X),1,sum)==5)
    cat(sprintf('Rows with no observed value: row'),i_allNA,'. Omitted\n')
    X = X[-i_allNA,]
    Y = Y[-i_allNA]
  }
  n=length(Y)


  X= as.matrix(cbind.data.frame(rep(1,nrow(X)),X))
  X.cc = as.matrix(na.omit(X))
  #indic=1:(ncol(X)-1)
  #Xnoint = X[,-1]
  #if(sum(indic)==0){rindic=0} else {rindic = as.matrix(is.na(X[,-1]))}
  #(TRUE=1 if missing, FALSE=0 if observed)
  rindic = as.matrix(is.na(X[,-1]))

  if(sum(rindic)>0){
    whichcolmissing = (1:ncol(rindic))[apply(rindic,2,sum)>0] #1 2 3
    whichcolXmissing = whichcolmissing+1 #2 3 4
    missingcols = length(whichcolmissing) #3
  }
  if(sum(rindic)==0){missingcols=0}

  if(missingcols>0){

    indicmiss = apply(as.matrix(rindic[,whichcolmissing]),1,sum)>0
    # 1=missing,0=observed for at least one of the variables for one observation

    Xnull = matrix(sapply(X, function(x) {x[is.na(x)] <- 0 ; x}),ncol=ncol(X))# 0<-NA
    rvx = apply(as.matrix(na.omit(X[,whichcolXmissing])),2,mean) #mean of the missing column()


    beta.start=beta.cc = glm(Y~X-1,family=binomial(link='logit'))$coef
    mu.start = mu.cc = apply(X.cc[,-1],2,mean)
    sig2.start = sig2.cc = var(X.cc[,-1])*(n-1)/n

    #sequence of beta
    seqbeta = matrix(NA,nrow=ncol(X),ncol=(maxruns+1)) # 5*1001
    #beta.start=glm(Y~X-1,family=binomial)$coef

    #sequence of sigma
    #seqsig.11=seqsig.12=seqsig.21=seqsig.22=matrix(NA,nrow=1,ncol=(maxruns+1))
    #seqsig=matrix(NA,nrow=length(whichcolXmissing),ncol=(maxruns+1)*length(whichcolXmissing))

    #sequence of amu
    #seqamu1 = seqamu2 = matrix(NA,nrow=(ncol(X)-missingcols),ncol=(maxruns+1))
    #seqamu = matrix(NA,nrow=length(whichcolXmissing),ncol=(maxruns+1)*(ncol(X)-1))

    #amu.start = matrix(0,length(whichcolXmissing),p)
    #for (i in 1:length(whichcolXmissing)){
    #  amu.start[i,] = lm(X[,whichcolXmissing[i]]~X[,-whichcolXmissing[i]]-1)$coef #x1 = amu[1,1]+ amu[1,2]*x2+ amu[1,3]*x3+ amu[1,4]*x4
    #}

    #mu.start = apply(X.cc[,-1],2,mean)
    #sig2.start = var(X.cc[,-1])*(n-1)/n
    #sig2.start=var(X[,whichcolXmissing],na.rm=T)#  var(x1,x2,x3)
    #sig2inv.start=solve(sig2.start)

    isum = imac = 0
    bigcounter = 0
    cstop=0.1 #tolerance
    mc.size=10

    while ((cstop>tol_em)*(imac<maxruns)|(imac<NBITERINIT)){ #stopping criteria for EM iterations
      rx = matrix(rep(0,mc.size*n*missingcols),nrow=n) #contains mc samples
      # E-step
      #	  t1=Sys.time()
      imac = imac + 1
      for(i in 1:n){
        if(indicmiss[i]){#NA exists in obs i
          idum= as.numeric(rindic[i,whichcolmissing]) # ex: 1 1 0, the missing pattern of missing var in ith obs
          patn= as.numeric(rindic[i,]) # ex: 1 1 0 0, the missing pattern of x1-x4 ith obs
          for(j1 in 1:mc.size){
            #if(missingcols==2){
            ###############################################?????????????????
            #rx[i,c(j1,mc.size+j1)] = Gibbs(iseed,Y,Xnull,beta.start,amu1.start,amu2.start,sig2inv.start,idum,patn,i,rvx)
            rx[i,seq(from=j1,to=(length(whichcolXmissing)-1)*mc.size+j1,by=mc.size)] = accept_reject(Y,Xnull,beta.start,mu.start,sig2.start,idum,patn,i,rvx)
            #if(missingcols==1){
            #  rx[i,j1] = Gibbs.1var(iseed,Y,Xnull,beta.start,amu1.start,sig2inv.start,idum,i,rvx,whichcolXmissing) }
          }
        } # end for(i...)
      }

      # M-step
      seqbeta[,imac] = beta.start
      startbeta = beta.start
      beta.optim = nlm(p=startbeta,f=fn.pvar,Y=Y,Xnull=Xnull,rx=rx,rindic=rindic,mc.size=mc.size,iterlim=20000, whichcolXmissing=whichcolXmissing)
      beta.start = beta.optim$estimate


      #inpute the missing value with the mean value of samples
      meanrx <- apply(rx[,1:mc.size],1,mean)
      for(i in 2:length(whichcolXmissing)){
        meanrx <- cbind(meanrx,apply(rx[,(mc.size*(i-1)+1):(mc.size*i)],1,mean))
        #meanrx<- cbind(apply(rx[,1:mrep],1,mean),apply(rx[,(mrep+1):(2*mrep)],1,mean))
      }
      newxvalues = Xnull
      newxvalues[,whichcolXmissing] = Xnull[,whichcolXmissing]+as.matrix(meanrx)*rindic[,whichcolmissing]

      #resid = matrix(0,length(whichcolXmissing),n)
      #for (i in 1:length(whichcolXmissing)){
      #  #amu.start = lm(newxvalues[,whichcolXmissing[1]]~newxvalues[,-whichcolXmissing]-1)$coef
      #  #res1 = lm(newxvalues[,whichcolXmissing[1]]~newxvalues[,-whichcolXmissing]-1)$resid
      #  amu.start[i,] = lm(newxvalues[,whichcolXmissing[i]]~newxvalues[,-whichcolXmissing[i]]-1)$coef #x1 = amu[1,1]+ amu[1,2]*x2+ amu[1,3]*x3+ amu[1,4]*x4
      #  resid[i,] = lm(newxvalues[,whichcolXmissing[i]]~newxvalues[,-whichcolXmissing[i]]-1)$resid
      #}
      #sig2.start = (resid%*%t(resid))/n #3x3

      mu.start = apply(newxvalues[,-1],2,mean)
      sig2.start = var(newxvalues[,-1])*(n-1)/n

      if (imac > 7) { #nb. iterations
        cstop=sum((beta.start-seqbeta[,imac-1])^2)
        if(print_iter==TRUE){
          print(cstop)
        }
      }

      if(imac>NBITERINIT){
        mc.size=mc.sizefin
      }
      if(print_iter==TRUE){
        cat(sprintf('iteration = %i ', imac))
        cat(sprintf('beta ='),beta.start,'\n')
        }
    }#end while ((cstop>0.0001)*(imac<maxruns)){ #stopping criteria for EM iterations
    #return(list(beta_est=beta_star))
    beta.est=seqbeta[,imac]

  }#end if(missingcols>0)

  if(missingcols==0){
    beta.est = glm(Y~X-1,family=binomial(link='logit'))$coef
    mu.start = apply(X[,-1],2,mean)
    sig2.start = var(X[,-1])*(n-1)/n
  }

  return(list(beta.est=beta.est,mu.est=mu.start,sig2.est=sig2.start))
}



#--------------------------generate data-----------------------------
generate.data = function(sample.size,beta.true,mu.true,sig2.true,perc.missing,iseed){
  #perc.missing : vector of missingness for each column
  set.seed(iseed)
  p = length(beta.true)-1
  X.orig = mvrnorm(sample.size,mu.true,sig2.true)
  X.orig = cbind(rep(1,sample.size),X.orig)

  pvector = exp(X.orig%*%(beta.true))/(1+exp(X.orig%*%beta.true))
  Y = rbinom(n=sample.size,size=1,prob=pvector)

  mis = matrix(0,p,sample.size)
  X=X.orig
  for(i in 1:p){
    mis[i,] = rbinom(n=sample.size,size=1,prob=1-perc.missing[i])
    mis[i,][mis[i,]==0]=NA
    X[,i+1] = X.orig[,i+1]*mis[i,]
  }

  datamatrix<- cbind(Y,X)
  data<- na.omit(datamatrix)

  return(list(X.orig=X.orig,Y.orig=Y,X=X,Y=Y,Y.cc=data[,1],X.cc=data[,-1]))
}

#--------------------------generate patterns-----------------------------
patterns = function(n){
  comb = NULL
  if (n<15) {
    for( i in 1:n) comb = rbind(cbind(1,comb),cbind(0,comb))
    return(comb)
  }
  else {error("this value will probably block your computer, try on your own risk")}
}


#--------------------------function to maximize in M-step-----------------------------

fn.pvar = function(beta.start,Y,Xnull,rx,rindic,mc.size,whichcolXmissing){
  # multinomial log likelihood
  # missingness can exist for all variables
  # first column of x consists of 1 only (intercept)
  # change newxvalues and columns in rindic to adapt to full missingness
  # rx = matrix(rep(0,mc.size*n*missingcols),nrow=n) #50*(500*3) contains mc samples
  # beta.optim = nlm(p=startbeta,f=fn.pvar,Y=Y,Xnull=Xnull,rx=rx,rindic=rindic,mc.size=mc.size,iterlim=20000, whichcolXmissing=whichcolXmissing)

  newxvalues = Xnull
  #mrep = ncol(rx)/2
  sumvec = rep(NA,mc.size)
  if(sum(as.numeric(rindic))>0){
    for(j in 1:mc.size) {
      newxvalues[,whichcolXmissing] = Xnull[,whichcolXmissing]+rx[,seq(from=j,to=(length(whichcolXmissing)-1)*mc.size+j,by=mc.size)]*rindic[,(whichcolXmissing-1)]
      #j:mc_size:length(whichcolXmissing)*mc.size+j
      #seq(from=j,to=length(whichcolXmissing)*mc.size+j,by=mc.size)

      ein = newxvalues%*%beta.start
      sumvec[j] = sum(ein*Y -log(1+exp(ein)))
    }
  }
  if(sum(as.numeric(rindic))==0){
    ein = xnull%*%beta.start
    sumvec = sum(ein*Y - log(1+exp(ein)))
  }
  return(Q1.beta.start = -mean(sumvec))
}


#--------------------------Accept-reject sampling (multi dim)-----------------------------
accept_reject=function(Y,Xnull,beta ,mu,sig2,idum,patn,i,rvx){
  # idum :missing pattern of missing variables 0 1 0
  # patn :missing pattern of x1-x4 0 1 0 0
  # rx = matrix(rep(0,mc.size*n*missingcols),nrow=n) #50*(500*3) contains mc samples
  # rvx = apply(as.matrix(na.omit(X[,whichcolXmissing])),2,mean) #mean of the missing column()
  # length(whichcolXmissing)=3


  #rx[i,seq(from=j1,to=(length(whichcolXmissing)-1)*mc.size+j1,by=mc.size)] = accept_reject(Y,Xnull,beta.start,mu.start,sig2.start,idum,patn,i,rvx)
  miss_col = which(patn==1)
  nb_miss_all = length(idum)
  output = rep(0,nb_miss_all) #3


  xi = Xnull[i,-1]
  x2 = xi[-miss_col]
  y = Y[i]

  mu1 = mu[miss_col]
  mu2 = mu[-miss_col]
  sigma11 = sig2[miss_col,miss_col]
  sigma12 = sig2[miss_col,-miss_col]
  sigma22 = sig2[-miss_col,-miss_col]
  sigma21 = sig2[-miss_col,miss_col]
  mu_cond = mu1+sigma12 %*% solve(sigma22)%*%(x2-mu2)
  sigma_cond = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21

  M=1
  u<-runif(1)*M
  z<-mvrnorm(n = 1, mu_cond, sigma_cond, tol = 1e-6, empirical = FALSE, EISPACK = FALSE) # draw from g
  #f<-fct_f(z,x1,y,beta,mu,sigma)
  #g<-dnorm(z,mu_2,sqrt(sigma_2),FALSE)
  #g<-mvn(z,mu_2,sigma_2,FALSE)
  while (u>fct_f(z,x2,y,beta ,mu,sig2,miss_col)/dmvnorm(z,mu_cond,sigma_cond,FALSE)/M){
    u<-runif(1)*M
    z<-mvrnorm(n = 1, mu_cond, sigma_cond, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  }
  #a single simulation of x2
  output[miss_col]<-z
  return(output)
  # rx[i,seq(from=j1,to=(length(whichcolXmissing)-1)*mc.size+j1,by=mc.size)] = accept_reject(Y,Xnull,beta.start,amu1.start,amu2.start,sig2inv.start,idum,patn,i,rvx)

}


fct_f <- function(x1,x2,y,beta,mu,sig2,miss_col){
  mu1 = mu[miss_col]
  mu2 = mu[-miss_col]
  sigma11 = sig2[miss_col,miss_col]
  sigma12 = sig2[miss_col,-miss_col]
  sigma22 = sig2[-miss_col,-miss_col]
  sigma21 = sig2[-miss_col,miss_col]
  x = rep(0,length(mu))
  x[miss_col]<-x1
  x[-miss_col]<-x2

  mu_cond = mu1+sigma12 %*% solve(sigma22)%*%(x2-mu2)
  sigma_cond = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21

  proba_mvn <- dmvnorm(x1, mean = mu_cond, sigma = sigma_cond, log = FALSE)
  proba_lr <- log_reg(y,c(1,x),beta,iflog=FALSE)
  return(proba_mvn * proba_lr)##JJ ici on a que le numerateur de f
}

log_reg <- function(y,x,beta,iflog=TRUE){
  res <- y*(x%*%beta) - log(1+exp(x%*%beta))
  if(iflog==TRUE)
    return(res)
  else
    return(exp(res))
}

