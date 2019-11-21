# This file contains the implementation of dti.smooth() for
# "dtiData" Adaptive smoothing in SE(3)

dwi.smooth <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dwi.smooth", function(object, ...) standardGeneric("dwi.smooth"))

setMethod("dwi.smooth", "dtiData", function(object,kstar,lambda=20,mask=NULL,kappa0=NULL,
                                            ncoils=1,sigma=NULL,level=NULL,
                                            vred=4,verbose=FALSE,dist=1,model=
                                              c("Gapprox","Gapprox2","Chi","Chi2")){
  ## check model
  model <- match.arg(model)
  args <- sys.call(-1)
  args <- c(object@call,args)
  sdcoef <- object@sdcoef
  level <- object@level
  vext <- object@voxelext[2:3]/object@voxelext[1]
  if(is.null(mask)) mask <- getmask(object,level)$mask
  nvoxel <- sum(mask)
  if(length(sigma)==1) {
    cat("using supplied sigma",sigma,"\n")
  } else {
    sigma <- numeric(object@ngrad)
    for(i in 1:object@ngrad){
      sigma[i] <- awssigmc(object@si[,,,i],12,mask,ncoils,vext,h0=1.25,verbose=verbose)$sigma
      cat("image ",i," estimated sigma",sigma[i],"\n")
    }
    cat("quantiles of estimated sigma values",quantile(sigma),"\n")
    sigma <- median(sigma)
    cat("using median estimated sigma",sigma,"\n")
  }
  model <- switch(model,"Chi2"=1,"Chi"=0,"Gapprox2"=2,"Gapprox"=3,1)
  #
  #  Chi2 uses approx of noncentral Chi^2 by rescaled central Chi^2 with adjusted DF
  #        and smoothes Y^2
  #  Chi uses approx of noncentral Chi^2 by rescaled central Chi^2 with adjusted DF
  #        and smoothes Y
  #   uses approximation of noncentral Chi by a Gaussian (with correct mean and variance)
  #        and smoothes Y
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad <- ngrad - ns0
  grad <- object@gradient[,-s0ind]
  bvalues <- object@bvalue[-s0ind]
  sb <- object@si[,,,-s0ind]
  s0 <- object@si[,,,s0ind]
  #
  #  rescale so that we have Chi-distributed values
  #
  if(ns0>1){
    dim(s0) <- c(prod(ddim),ns0)
    if(model%in%c(1,2)){
      s0 <- sqrt(s0^2%*%rep(1,ns0))
    } else {
      s0 <- s0%*%rep(1,ns0)
    }
  }
  #
  #  keep only voxel in mask
  #
  s0 <- s0[mask]
  dim(sb) <- c(prod(ddim),ngrad)
  sb <- sb[mask,]
    # S0 will be noncentral Chi with 2*ns0*ncoils DF for model 1 and 2
  sb <- sb/sigma
  s0 <- s0/sigma
  if(model==1)    sigma <- sigma^2

  z <- aws::smse3(sb, s0, bvalues, grad, ns0, kstar, lambda,
                  kappa0, mask, vext, vred, ncoils, model, dist, verbose=verbose)

  ngrad <- ngrad+1
  si <- array(0,c(prod(ddim),ngrad))
  si[mask,1] <- z$th0
  si[mask,-1] <- z$th
  dim(si) <- c(ddim,ngrad)
  #
  #  back to original scale
  #
  s0factor <- switch(model+1,ns0,ns0,sqrt(ns0),ns0)
  bvalue <- c(0,object@bvalue[-object@s0ind])
  si[,,,1] <-  th0*sigma/s0factor
  si[,,,-1] <- z$th*sigma
  object@si <- if(model==1) sqrt(si) else si
  object@gradient <- grad <- cbind(c(0,0,0),grad)
  object@bvalue <- bvalue
  object@btb <- sweep(create.designmatrix.dti(grad), 2, bvalue, "*")
  object@s0ind <- as.integer(1)
  object@replind <- as.integer(1:ngrad)
  object@ngrad <- as.integer(ngrad)
  object@call <- args
  object
}
)
