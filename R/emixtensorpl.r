#
#  spatially extended tensor mixture models
#
#dwiMixtensor <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

#setGeneric("dwiMixtensor", function(object,  ...) standardGeneric("dwiMixtensor"))

#setMethod("dwiMixtensor","dtiData",
dwiExtMt <- function(object, ncomp=2, 
               neighbors=matrix(c(0,0,0,1,-1,0,0,.2,1,0,0,.2,
                                  0,-1,0,.2,0,1,0,.2,0,0,-1,.2,
                                  0,0,1,.2),4,7), 
          model=c("MTiso"), 
          reltol=1e-5, maxit=500,ngc=1000, nguess=100*ncomp^2, 
          msc=c("BIC","AIC","AICC","none"), mc.cores = setCores(,reprt=FALSE)){
#
#  row's of matrix neighbors give relative x,y,z coordinates plus weights 
#  ncomp refers to the number of compartments, no model comparisons here
#        different models need to be evaluated seperately
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
#  Optimization methods: 
#     L-BFGS-B for tensor mixture models with isotropic compartment
#
  ## check model
  model <- match.arg(model)
  ## check msc
  msc <- match.arg(msc)
  factr <- reltol/1e-14 ## this is 1e6 
  set.seed(1)
  bvalue <- object@bvalue[-object@s0ind]
  maxbv <- max(bvalue)
  bvalue <- bvalue/maxbv
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  nvox <- prod(ddim)
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(5*(1+3*ncomp)>ngrad0){
#     maxcomp <- max(1,trunc((ngrad0-5)/15))
     cat("number of components reduced to", ncomp,"due to insufficient
          number of gradient directions\n")
  }
#
#   which model should be used
#
  imodel <- 2
  alpha <- 0
  lambda <- 0 
#
#  First tensor estimates to generate eigenvalues and -vectors
#
  prta <- Sys.time()
  cat("Start tensor estimation at",format(prta),"\n")
  tensorobj <- dtiTensor(object, mc.cores = mc.cores)
  cat("Start evaluation of eigenstructure at",format(Sys.time()),"\n")
  z <- dtieigen(tensorobj@D, tensorobj@mask, mc.cores = mc.cores)
  rm(tensorobj)
  gc()
  fa <- array(z$fa,ddim)
  ev <- array(z$ev,c(3,ddim))*maxbv
#
#  rescale by bvalue to go to implemented scale
#
  andir <- array(z$andir,c(3,2,ddim))
  rm(z)
  gc()
#
#  prepare parameters for searching initial estimates
#
   lambdahat <- ev[3,,,] 
# use third ev instead of (ev[2,,,]+ev[3,,,])/2 to avoid effects from mixtures
   alphahat <- if(imodel==2) (ev[1,,,]-lambdahat)/lambdahat else alpha
   lambdahat <- median(lambdahat[!is.na(alphahat)&fa>.3])
   if(imodel==2) alphahat <- median(alphahat[!is.na(alphahat)&fa>.3])
   fahat <- alphahat/sqrt(3+2*alphahat+alphahat^2)
   cat("Using lambda_2=",lambdahat,"fa=",fahat," and alpha=",alphahat,"in initial estimates\n")
  cat("Start search outlier detection at",format(Sys.time()),"\n")
#
#  replace physically meaningless S_i by mena S_0 values
#
  z <- sioutlier(object@si,s0ind,mc.cores=mc.cores)
#
#  this does not scale well with openMP
#
  cat("End search outlier detection at",format(Sys.time()),"\n")
  si <- array(z$si,c(ngrad,ddim))
  index <- z$index
  rm(z)
  gc()
  cat("Start generating auxiliary objects",format(Sys.time()),"\n")
#
#  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
#
  nvox <- prod(ddim[1:3])
  cat("sweeps0:")
  t1 <- Sys.time()
  if(mc.cores==1||ngrad0>250){
  z <- .Fortran("sweeps0",# mixtens.f
                as.double(si[-s0ind,,,,drop=FALSE]),
                as.double(si[s0ind,,,,drop=FALSE]),
                as.integer(nvox),
                as.integer(ns0),
                as.integer(ngrad0),
                as.integer(object@level),
                siq=double(nvox*ngrad0),
                s0=double(nvox),
                vsi=double(nvox),
                mask=logical(nvox),
                DUPL=FALSE,
                PACKAGE="dti")[c("siq","s0","vsi","mask")]
  t2 <- Sys.time()
  cat(difftime(t2,t1),"for",nvox,"voxel\n")
  s0 <- array(z$s0,ddim[1:3])
  siq <- array(z$siq,c(ngrad0,ddim[1:3]))
#
#  siq is permutated c(4,1:3)
#
  sigma2 <- array(z$vsi,ddim[1:3])
  mask <- array(z$mask,ddim[1:3])
  } else {
  mc.cores.old <- setCores(,reprt=FALSE)
  setCores(mc.cores)
  z <- matrix(.Fortran("sweeps0p",# mixtens.f
                as.double(si[-s0ind,,,,drop=FALSE]),
                as.double(si[s0ind,,,,drop=FALSE]),
                as.integer(nvox),
                as.integer(ns0),
                as.integer(ngrad0),
                as.integer(object@level),
                siq=double(nvox*(ngrad0+3)),
                as.integer(ngrad0+3),
                DUPL=FALSE,
                PACKAGE="dti")$siq,ngrad0+3,nvox)
  t2 <- Sys.time()
  cat(difftime(t2,t1),"for",nvox,"voxel\n")
  setCores(mc.cores.old,reprt=FALSE)
  s0 <- array(z[ngrad0+1,],ddim[1:3])
  siq <- array(z[1:ngrad0,],c(ngrad0,ddim[1:3]))
#
#  siq is permutated c(4,1:3)
#
  sigma2 <- array(z[ngrad0+2,],ddim[1:3])
  mask <- array(as.logical(z[ngrad0+3,]),ddim[1:3])
  }
  rm(si)
  rm(z)
  gc()
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  nn <- dim(neighbors)[2] # number of voxel in neighborhood
  z <- .Fortran("nvindex",
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.logical(mask),
                as.integer(dim(neighbors)[2]),
                as.integer(neighbors[-4,]),
                ind=integer(nvox*nn),
                nind=as.integer(nvox),
                DUPL=TRUE,PACKAGE="dti")[c("ind","nind")]
  nmask <- z$nind
  ind <- matrix(z$ind,nn,nmask)
#
# ind contains the one-dimensional index of voxel from the neighborhood 
# for each voxel in mask. Index 0 corresponds to voxel outside the cube or mask
#
  npar <- 2+(2+nn)*ncomp
 #
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
#
#   determine initial estimates for orientations 
#
  cat("Start search for initial directions at",format(Sys.time()),"\n")
  data("polyeders", envir = environment())
  polyeder <- icosa3
  vert <- polyeder$vertices
# remove redundant directions
  vind <- rep(TRUE,dim(vert)[2])
  vind[vert[1,]<0] <- FALSE
  vind[vert[1,]==0 & vert[2,] <0] <- FALSE
  vind[vert[1,]==0 & vert[2,] == 0 &vert[3,]<0] <- FALSE
  vert <- vert[,vind]
  cat("End generating auxiliary objects",format(Sys.time()),"\n")
#
#  compute initial estimates (EV from grid and orientations from icosa3$vertices)
#
  dim(siq) <- c(ngrad0,nvox)

  siind <- matrix(0,ncomp+1,nvox)
  siind[1,!mask] <- -1
  krit <- numeric(nvox)
  nvoxm <- sum(mask)
  if(mc.cores<=1){
    z <-   getsiind2(siq[,mask],sigma2[mask],grad,bvalue,t(vert),alphahat,lambdahat,
                       ncomp,maxc=maxc,nguess=nguess)
    krit[mask] <- z$krit # sqrt(sum of squared residuals) for initial estimates
    siind[,mask] <- z$siind # components 1: model order 2: 
                       # grid index for EV 2+(1:m) index of orientations
  } else {
    mc.cores.old <- setCores(,reprt=FALSE)
    setCores(mc.cores)
    x <- array(0,c(ngrad0+1,nvoxm))
    x[1:ngrad0,] <- siq[,mask]
    x[ngrad0+1,] <- sigma2[mask]
    nvico <- dim(vert)[2]
    dgrad <- matrix(abs(grad%*%vert),ngrad0,nvico)
    dgrad <- dgrad/max(dgrad)
    dgradi <- matrix(abs(t(vert)%*%vert),nvico,nvico)
    dgradi <- dgradi/max(dgradi)
    isample <- selisample(nvico,ncomp,nguess,dgradi,maxc)
    nguess <- length(isample)/ncomp
    cat("using ",nguess,"guesses for initial estimates\n")
    z <- plmatrix(x,pgetsiind2,grad=grad,bv=bvalue,nvico=nvico,
                  dgrad=dgrad,dgradi=dgradi,isample=isample,alpha=alphahat,
                  lambda=lambdahat,maxcomp=ncomp,maxc=maxc,nguess=nguess)
    setCores(mc.cores.old,reprt=FALSE)
    krit[mask] <- z[1,] # risk
    siind[,mask] <- z[-1,] # siind 
  }
  cat("Model orders for initial estimates")
  print(table(siind[1,]))
  cat("End search for initial values at",format(Sys.time()),"\n")
#  logarithmic eigen values
  orient <- array(0,c(2,ncomp,ddim))
  prt0 <- Sys.time()
#
#   loop over voxel in volume
#
  siind <- siind[-1,]
  vertorient <- array(.Fortran("parofor",
                      as.double(vert),
                      as.integer(dim(vert)[2]),
                      orient=double(2*dim(vert)[2]),
                      DUPL=FALSE,
                      PACKAGE="dti")$orient,c(2,dim(vert)[2]))
  lambdaest <- alphaest <- sigma2 <- numeric(nvox)
  orient <- array(0,c(2,ncomp,nvox))
  mix <- matrix(0,ncomp,nvox)
  lower <- c(0,0.6,rep(c(1e-3,-Inf,-Inf),ncomp),rep(1e-3,ncomp*(nn-1)))
  upper <- c(1e2/max(bvalue),10,rep(c(1e3,Inf,Inf),ncomp),rep(1e3,ncomp*(nn-1)))
  cat("Start optimization",format(Sys.time()),"\n")
  
  for (i in 1:nmask){
      param <- c(lambdahat,alphahat)
      for(j in 1:ncomp){
         param <- c(param,1,vertorient[,siind[j,ind[1,i]]])
      }
      indi <- ind[,i]
      indpos <- indi>0
      indi <- indi[indpos]
      param <- c(param,rep(1,ncomp*(length(indi)-1)))
      npar <- length(param)
      siqi <- siq[,indi]
      z <- optim(param,erskmxl2,edrskml2,siq=siqi,grd=t(grad),b=bvalue,vw=neighbors[4,indpos],
                 method="L-BFGS-B",lower=lower[1:npar],
                 upper=upper[1:npar],control=list(maxit=maxit,reltol=reltol))
      param <- z$par
      sigma2[indi[1]] <- z$value
#      cat("value",z$value,"param",param,"\n")
      lambdaest[indi[1]] <- param[1]
      alphaest[indi[1]] <- param[2]
      orient[,,indi[1]] <- param[rep(1:ncomp,rep(2,ncomp))*3+1:2]
      mix[,indi[1]] <- param[3*(1:ncomp)]/(1+sum(param[3*(1:ncomp)]))
      if(i%/%100*100==i) cat("processed",i,"out of",nmask,"voxel --> time:",
                               format(Sys.time()),"\n")
  }
  cat("End optimization",format(Sys.time()),"\n")
  dim(sigma2) <- ddim
  order <- array(0, ddim)
  order[mask] <- ncomp
  lev <- matrix(0,2,nvox)
  lev[2,] <- lambdaest
  lev[1,mask] <- alphaest*lambdaest
  dim(lev) <- c(2,ddim)
  dim(mix) <- c(ncomp, ddim)
  dim(orient) <- c(2,ncomp, ddim)
  #
  method <- "MTiso"
  model <- "iso-prolate"
  invisible(new("dwiMixtensor",
                model = model,
                call   = args,
                ev     = lev/maxbv,
                mix    = mix,
                orient = orient,
                order  = order,
                p      = 0,
                th0    = s0,
                sigma  = sigma2,
                scorr  = array(1,c(1,1,1)), 
                bw     = c(0,0,0), 
                mask   = mask,
                hmax   = 1,
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb    = object@btb,
                ngrad  = object@ngrad, # = dim(btb)[2]
                s0ind  = object@s0ind,
                replind = object@replind,
                ddim   = object@ddim,
                ddim0  = object@ddim0,
                xind   = object@xind,
                yind   = object@yind,
                zind   = object@zind,
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                rotation = object@rotation,
                source = object@source,
                outlier = index,
                scale = 1,
                method = method)
            )
}
#
#   compution the risk
#

erskmxl2 <- function(par,siq,grd,b,vw){
   npar <- length(par)
   ng <- dim(grd)[2]
   nvox <- length(vw)
   risk <- .Fortran("erskmxl2",
                 as.double(par),
                 as.integer(npar),
                 as.double(siq),
                 as.double(grd),
                 as.double(b),
                 as.integer(ng),
                 double(nvox),
                 as.double(vw),
                 as.integer(nvox),
                 risk=double(1),
                 DUPL=FALSE,
                 PACKAGE="dti")$risk
#   cat("risk",risk,"\n")
   risk
}
#
#   compute gradient with respect to parameters
#
edrskml2 <- function(par,siq,grd,b,vw){
  npar <- length(par)
  ng <- dim(grd)[2]
  nvox <- length(vw)
  ncomp <- (npar-2)/(nvox+2)
  gradient <- .Fortran("edrskml2",
                   as.double(par),
                   as.integer(npar),
                   as.integer(ncomp),
                   as.double(siq),
                   as.double(grd),
                   as.double(b),
                   as.integer(ng),
                   double(nvox),#fval
                   as.double(vw),
                   as.integer(nvox),
                   double(3*ncomp*nvox),#dval
                   double(3*ncomp*nvox),#drisk0
                   double(2*nvox),#driskla
                   gradient=double(npar),# drisk
                   double(nvox),#dlam
                   double(nvox),#dalpha
                   DUPL=FALSE,
                   PACKAGE="dti")$gradient
  gradient
}
