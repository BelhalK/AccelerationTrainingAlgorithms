#' @param model an object of class SaemixModel, created by a call to the
#' function \code{\link{saemixModel}}
#' @param data an object of class SaemixData, created by a call to the function
#' \code{\link{saemixData}}
#' @param control a list of options, see \code{\link{saemixControl}}
#' @return An object of class SaemixObject containing the results of the fit of
#' the data by the non-linear mixed effect model. A summary of the results is
#' printed out to the terminal, and, provided the appropriate options have not
#' been changed, numerical and graphical outputs are saved in a directory.
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{SaemixObject}}, \code{\link{saemixControl}},
#' \code{\link{plot.saemix}}
#' @references Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear
#' mixed effects models. Computational Statistics and Data Analysis 49, 4
#' (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
#' 20th meeting of the Population Approach Group in Europe, Athens, Greece
#' (2011), Abstr 2173.
#' @keywords models
#' @examples
#' 
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
#'    name.group=c("Id"),name.predictors=c("Dose","Time"),
#'    name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'    units=list(x="hr",y="mg/L", covariates=c("kg","-")), name.X="Time")
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' 
#' # Not run (strict time constraints for CRAN)
#' # saemix.fit<-saemix(saemix.model,saemix.data,list(seed=632545,directory="newtheo",
#' # save=FALSE,save.graphs=FALSE))
#' 
#' # Prints a summary of the results
#' # print(saemix.fit)
#' 
#' # Outputs the estimates of individual parameters
#' # psi(saemix.fit)
#' 
#' # Shows some diagnostic plots to evaluate the fit
#' # plot(saemix.fit)
#' 
#' 
#' @export saemix
saemix_post_cat<-function(model,data,control=list()) {

 if(class(model)!="SaemixModel") {
    cat("Please provide a valid model object (see the help page for SaemixModel)\n")
    return()
  }
  if(class(data)!="SaemixData") {
    cat("Please provide a valid data object (see the help page for SaemixData)\n")
    return()
  }
  

  saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
#  saemixObject<-new(Class="SaemixObject",data=saemix.data, model=saemix.model,options=saemix.options)
  opt.warn<-getOption("warn")
  if(!saemixObject["options"]$warnings) options(warn=-1)

  saemix.options<-saemixObject["options"]
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
  saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
  saemix.data@ntot.obs<-dim(saemix.data@data)[1]
#  showall(saemixObject)

# Initialising random generator
  OLDRAND<-TRUE
  set.seed(saemix.options$seed)

############################################
#  Main Algorithm
############################################
  
# Initialisation - creating several lists with necessary information extracted (Uargs, Dargs, opt,varList, suffStat)
  xinit<-initialiseMainAlgo_cat(saemix.data,saemix.model,saemix.options)
  
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
  theta0<-c(fixed.psi,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res])

  parpop<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1),ncol=(Uargs$nb.parameters+length(Uargs$i1.omega2)))
  colnames(parpop)<-c(saemix.model["name.modpar"], saemix.model["name.random"])
  allpar<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1), ncol=(Uargs$nb.betas+length(Uargs$i1.omega2)))
  colnames(allpar)<-c(saemix.model["name.fixed"],saemix.model["name.random"])
  
  parpop[1,]<-theta0
  allpar[1,]<-xinit$allpar0
  
  # using several Markov chains - only useful if passed back to main routine...
  #   chdat<-new(Class="SaemixRepData",data=saemix.data, nb.chains=saemix.options$nb.chains)
  #   NM<-chdat["NM"]
  #   IdM<-chdat["dataM"]$IdM
  #   yM<-chdat["dataM"]$yM
  #   XM<-chdat["dataM"][,saemix.data["name.predictors"],drop=FALSE]
  
# List of sufficient statistics - change during call to stochasticApprox
  suffStat<-list(statphi1=0,statphi2=0,statphi3=0,statrese=0)
  phi<-array(data=0,dim=c(Dargs$N, Uargs$nb.parameters, saemix.options$nb.chains))

# structural model, check nb of parameters
  structural.model<-saemix.model["model"]
  #  nb.parameters<-saemix.model["nb.parameters"]

  xmcmc<-estep_cat(1, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM, saemixObject)

  # xmcmc<-estep_newkernel(1, Uargs, Dargs, opt, structural.model, mean.phi, varList, DYF, phiM)
  # varList<-xmcmc$varList
  DYF<-xmcmc$DYF
  phiM<-xmcmc$phiM
  post_rwm<-xmcmc$post

  post_new<-xmcmc$post_new

  return(list(post_rwm = post_rwm,post_new = post_new))

}
