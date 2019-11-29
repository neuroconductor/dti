#
#
#      estimate variance parameter in a multicoil system
#
#
awslsigmc <- function(y,                 # data
                      steps,             # number of iteration steps for PS
                      mask = NULL,       # data mask, where to do estimation
                      ncoils = 1,        # number of coils for parallel MR image acquisition
                      vext = c( 1, 1),   # voxel extensions
                      lambda = 5,       # adaptation parameter for PS
                      minni = 2,         # minimum sum of weights for estimating local sigma
                      hsig = 5,          # bandwidth for median smoothing local sigma estimates
                      sigma = NULL,
                      family = c("NCchi"),
                      verbose = FALSE,
                      trace = FALSE,
                      u=NULL#,
                      #bc=FALSE # bias correction ...
) {
## Code has been moved to aws package
## function from aws is called for compatibility reasons
aws::awsLocalSigma(y, steps, mask, ncoils, vext,
      lambda, minni, hsig, sigma, family, verbose, trace, u)
}


#
#
#      estimate variance parameter in a multicoil system
#
#
awssigmc <- function(y,                 # data
                     steps,             # number of iteration steps for PS
                     mask = NULL,       # data mask, where to do estimation
                     ncoils = 1,        # number of coils for parallel MR image acquisition
                     vext = c( 1, 1),   # voxel extensions
                     lambda = 20,       # adaptation parameter for PS
                     h0 = 2,            # initial bandwidth for first step in PS
                     verbose = FALSE,
                     sequence = FALSE,  # return estimated sigma for intermediate steps of PS?
                     hadj = 1,          # adjust parameter for density() call for mode estimation
                     q = .25,  # for IQR
                     qni = .8,
                     method=c("VAR","MAD")  # for variance, alternative "MAD" for mean absolute deviation
) {
  method <- switch(method,"VAR"="awsVar","MAD"="awsMAD")
  aws::estGlobalSigma(y, mask, ncoils, steps, vext, lambda, h0,
                   hadj, q, qni, sequence=sequence, method=method)
}

aflsigmc <- function(y,ncoils,level=NULL,mask=NULL,h=2,hadj=1,vext = c( 1, 1)){
  ##
  ##   estimate effective sigma and effective ncoils (L) according to Aja-Fernandez 2013
  ##
  aws::AFLocalSigma(y,ncoils,level,mask,h,hadj,vext)
}

afsigmc <- function(y,                 # data
                    level = NULL,             # threshold for background separation
                    mask = NULL,       # data mask, where to do estimation, needs to refer to background if level == NULL
                    ncoils = 1,        # number of coils for parallel MR image acquisition
                    vext = c( 1, 1),   # voxel extensions
                    h = 2,             # bandwidth for local averaging
                    verbose = FALSE,
                    hadj = 1,           # adjust parameter for density() call for mode estimation
                    method = c("modevn","modem1chi","bkm2chi","bkm1chi")# methods according to table 2 in Aja-Ferbnandez (2009)
) {
  method <- switch(method, "modevn"="AFmodevb",
                           "modem1chi"="AFmodevb",
                           "bkm2chi"="AFbkm2chi",
                           "bkm1chi"="AFbkm1chi")
  aws::estGlobalSigma(y, mask, ncoils, vext=vext,
                      lambda=NULL, h0=h, hadj=hadj,
                      level=level, method="method")
}
