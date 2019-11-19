# This file contains the implementation of dti.smooth() for
# "dtiData" Adaptive smoothing in SE(3) considering b=0 as an individual shell

dwi.smooth.ms <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dwi.smooth.ms", function(object, ...) standardGeneric("dwi.smooth.ms"))

setMethod("dwi.smooth.ms",
          "dtiData",
          function(object,
                   kstar,
                   lambda = 12,
                   kappa0 = .5,
                   ncoils = 1,
                   sigma = NULL,
                   ws0 = 1,
                   level = NULL,
                   mask = NULL,
                   xind = NULL,
                   yind = NULL,
                   zind = NULL,
                   verbose = FALSE,
                   usemaxni = TRUE,
                   memrelease = TRUE){

            # make the call part of the object
            args <- sys.call(-1)
            args <- c(object@call,args)

            # we need a lot of object properties
            ddim <- object@ddim
            ngrad <- object@ngrad
            s0ind <- object@s0ind
            ns0 <- length(s0ind)
            ngrad <- ngrad - ns0
            grad <- object@gradient[,-s0ind]
            bvalues <- object@bvalue[-s0ind]
            sdcoef <- object@sdcoef
            level <- object@level
            vext <- object@voxelext[2:3]/object@voxelext[1]


            if (length(sigma) == 1) {
              sigmacase <- 1
              if (verbose) cat("using supplied sigma ", sigma, "\n")
            } else if (identical(dim(sigma), ddim)) {
              sigmacase <- 2
              if (verbose) cat("using supplied array of sigma\n")
            } else if (identical(dim(sigma), c(ddim, nshell+1L))) {
              sigmacase <- 3
              if (is.null(mask)) mask <- getmask(object, level)$mask
              if (verbose) cat("using supplied array of sigma per shell\n")
            } else {
              if (is.null(mask)) mask <- getmask(object, level)$mask
              sigma <- numeric(object@ngrad)
              for(i in 1:object@ngrad){
                ## FOR FUTURE USE:
                # awslsigmc(object@si[ , , , i], 12, ncoils = ncoils, lambda = 5, verbose = verbose, hsig = 5, mask = mask)
                ## END
                sigma[i] <- awssigmc(object@si[,,,i], 12, mask, ncoils, vext, h0=1.25, verbose=verbose)$sigma
                if (verbose) cat("image ", i," estimated sigma", sigma[i], "\n")
              }
              if (verbose) cat("quantiles of estimated sigma values", quantile(sigma), "\n")
              sigma <- median(sigma)
              if (verbose) cat("using median estimated sigma", sigma,"\n")
              sigmacase <- 1
            }


            #
            sb <- object@si[,,,-s0ind]
            s0 <- object@si[,,,s0ind]
            if(is.null(kappa0)){
              #  select kappa based on variance reduction on the sphere
              warning("You need to specify  kappa0  returning unsmoothed object")
              return(object)
            }
            #
            #  rescale so that we have Chi-distributed values
            #
            if (sigmacase == 1) {
              sb <- sb/sigma
              s0 <- s0/sigma
            } else if (sigmacase == 2) {
              sb <- sweep(sb, 1:3, sigma, "/")
              s0 <- sweep(s0, 1:3, sigma, "/")
            } else { # now is per shell
              s0 <- sweep(s0, 1:3, sigma[, , , 1], "/")
              for (shnr in 1:nshell) {
                indbv <- (1:length(bv))[bv == ubv[shnr]]
                sb[, , , indbv] <- sweep(sb[, , , indbv], 1:3, sigma[, , , shnr+1], "/")
              }
            }
            if(ns0>1){
              dim(s0) <- c(prod(ddim),ns0)
              s0 <- s0%*%rep(1/sqrt(ns0),ns0)
              #  make sure linear combination of s0 has same variance as original
              dim(s0) <- ddim
            }
            if (is.null(mask)) mask <- getmask(object, level)$mask
            # determine minimal subcube that contains the mask ore use supplied information
            ##
            ##  check memory size needed for largest vector
            ##
            vsize <- (nshell+1)*ngrad*sum(mask)
            if(vsize>2^31-1){
              stop("region specified is to large, a vector of size",vsize,"needs to be passed through
      the .Fortran, this is limited to 2^31-1 in R\n please use a smaller mask")
            }
            #
            #  reduce data to contain only voxel from brain mask
            #
            dim(sb) <- c(prod(ddim),ngrad)
            s0 <- s0[mask]
            sb <- sb[mask,]

            z <- aws::smse3ms(sb, s0, bv, grad, ns0, kstar, kappa0,
                            mask, vext=vext, ncoils=ncoils, level=level,
                            verbose=verbose, usemaxni=usemaxni)

            #
            #  one s0 image only
            #
            ngrad <- ngrad+1
            si <- array(0,c(prod(ddim),ngrad))
            si[mask,1] <- z$th0
            si[mask,-1] <- z$th
            dim(si) <- c(ddim,ngrad)
            #
            #  back to original scale
            #
            if (length(sigma) > 1) {
              if(length(dim(sigma)) == 3) {
                si[ , , , 1] <-  z$th0/sqrt(ns0)*sigma
              } else {
                si[, , , 1] <-  z$th0/sqrt(ns0)*sigma[, , , 1]
              }
            } else {
              si[, , , 1] <-  z$th0/sqrt(ns0)*sigma
            }
            #  go back to original s0 scale
            #  for the DWI we need to scale back differently
            if (sigmacase == 1) {
              si[, , , -1] <- z$th*sigma
            } else if (sigmacase == 2) {
              si[, , , -1] <- sweep(z$th, 1:3, sigma, "*")
            } else if (sigmacase == 3) { # now sigma is an array with (identical(dim(sigma), ddim) == TRUE)
              for (shnr in 1:nshell) {
                indth <- (1:length(bv))[bv == ubv[shnr]]
                indsi <- indth+1
                si[, , , indsi] <- sweep(z$th[, , , indth], 1:3, sigma[, , , shnr+1], "*")
                }
            } else {
                si[xind, yind, zind, -1] <- z$th*sigma[,,,-1]
            }
            object@si <-  si
            object@gradient <- grad <- cbind(c(0,0,0),grad)
            object@bvalue <- bvalue <- c(0,object@bvalue[-object@s0ind])
            object@btb <- sweep(create.designmatrix.dti(grad), 2, bvalue, "*")
            object@s0ind <- as.integer(1)
            object@replind <- as.integer(1:ngrad)
            object@ngrad <- as.integer(ngrad)
            object@call <- args
            object
          }
)
