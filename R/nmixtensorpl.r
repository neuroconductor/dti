dwiMixtensor2 <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor2", function(object,  ...) standardGeneric("dwiMixtensor2"))

setMethod("dwiMixtensor2","dtiData",function(object, maxcomp=3, fa=NULL, lambda=NULL, factr=1e7, maxit=5000,ngc=1000, nguess=100*maxcomp^2,msc="BIC",thinit=NULL, 
    mc.cores = setCores(,reprt=FALSE)){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
#  uses L-BFGS-B  
#  models with isotropic compartment
#     if(!is.null(FA)&!is.null(lambda)) use max. EV / FA as provided
#     if(!is.null(FA)&is.null(lambda)) use FA as provided
#     if(is.null(FA)) estimate both EV's
#  lambda corresponds to maximum EV
  set.seed(1)
#
#  pre-processing
#
  theta <- .5
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  nvox <- prod(ddim)
  s0ind <- object@s0ind
  bvalue <- object@bvalue[-s0ind]
  maxbv <- max(bvalue)
  bvalue <- bvalue/maxbv
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(5*(1+3*maxcomp)>ngrad0){
#     maxcomp <- max(1,trunc((ngrad0-5)/15))
     cat("Maximal number of components reduced to", maxcomp,"due to insufficient
          number of gradient directions\n")
  }
#
#   which model should be used
#
if(is.null(fa)){
   model <- 2  
} else {
   fa2 <- fa^2
   alpha <- (fa2+sqrt(fa2*(3-2*fa2)))/(1-fa2)
   if(is.null(lambda)){
      model <- 1
   } else {
      model <- 0
      lambda <- lambda*maxbv/(1+alpha)
   }
}
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
  lambdahat <- ev[3,,,] 
# use third ev instead of (ev[2,,,]+ev[3,,,])/2 to avoid effects from mixtures
  alphahat <- if(model==2) (ev[1,,,]-lambdahat)/lambdahat else alpha
  lambdahat <- median(lambdahat[!is.na(alphahat)&fa>0.3])
  if(model==2) alphahat <- median(alphahat[!is.na(alphahat)&fa>.3])
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
  npar <- model+3*(0:maxcomp)
#
#   compute penalty for model selection, default BIC
#
  penIC <- switch(msc,"AIC"=2*npar/ngrad0,"BIC"=log(ngrad0)*npar/ngrad0,
                  "AICC"=(1+npar/ngrad0)/(1-(npar+2)/ngrad0),
                  "None"=log(ngrad0)-log(ngrad0-npar),
                  log(ngrad0)*npar/ngrad0)
  cat("End generating auxiliary objects",format(Sys.time()),"\n")
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
#
#   determine initial estimates for orientations 
#
  cat("Start search for initial directions at",format(Sys.time()),"\n")
  data("polyeders")
  polyeder <- icosa3
  vert <- polyeder$vertices
# remove redundant directions
  vind <- rep(TRUE,dim(vert)[2])
  vind[vert[1,]<0] <- FALSE
  vind[vert[1,]==0 & vert[2,] <0] <- FALSE
  vind[vert[1,]==0 & vert[2,] == 0 &vert[3,]<0] <- FALSE
  vert <- vert[,vind]
#
#  compute initial estimates (EV from grid and orientations from icosa3$vertices)
#
  siind <-   getsiind2(siq,mask,sigma2,grad,bvalue,t(vert),alphahat,lambdahat,
                       maxcomp,maxc=maxc,nguess=nguess)
  krit <- siind$krit # sqrt(sum of squared residuals) for initial estimates
  siind <- siind$siind # components 1: model order 2: 
                       # grid index for EV 2+(1:m) index of orientations
  cat("Model orders for initial estimates")
  print(table(siind[1,,,]))
  cat("End search for initial values at",format(Sys.time()),"\n")
  siind <- siind[-1,,,]
# only need indices of initial directions
  orient <- array(0,c(2,maxcomp,ddim))
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  prt0 <- Sys.time()
#
#     C-Code
#
if(mc.cores<=1){
  cat("Starting parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
  dim(siq) <- c(ngrad0,nvox)
  dim(siind) <- c(maxcomp,nvox)
  nvoxm <- sum(mask)
  z <- switch(model+1,.C("mixtrl0", 
          as.integer(nvoxm),#n1
          as.integer(siind[,mask]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bvalue),#bv_in
          lambda  = as.double(lambda),#lambda_in
          alpha   = as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(sigma2[mask]),#sigma2
          as.double(vert),#vert
          as.double(siq[,mask]),#siq_in
          sigma2  = double(nvoxm),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvoxm),#order_ret selected order of mixture
          mix     = double(maxcomp*nvoxm),#mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")],
          .C("mixtrl1", 
          as.integer(nvoxm),#n1
          as.integer(siind[,mask]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bvalue),#bv_in
          as.double(lambdahat),#lambda_in
          alpha = as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(sigma2[mask]),#sigma2
          as.double(vert),#vert
          as.double(siq[,mask]),#siq_in
          sigma2  = double(nvoxm),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvoxm),#order_ret selected order of mixture
          lambda  = double(nvoxm),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvoxm),#mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")],
          .C("mixtrl2", 
          as.integer(nvoxm),#n1
          as.integer(siind[,mask]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bvalue),#bv_in
          as.double(lambdahat),#lambda_in
          as.double(alphahat),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(sigma2[mask]),#sigma2
          as.double(vert),#vert
          as.double(siq[,mask]),#siq_in
          sigma2  = double(nvoxm),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvoxm),#order_ret selected order of mixture
          alpha   = double(nvoxm),#alpha_ret alpha=(lambda_1-lambda_2)/lambda_2 
          lambda  = double(nvoxm),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvoxm),#mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")])
  cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
  sigma2 <-  array(0,ddim)
  sigma2[mask] <- z$sigma2
  orient <- matrix(0,2*maxcomp,nvox)
  orient[,mask] <- z$orient
  dim(orient) <- c(2, maxcomp, ddim)
  order <- array(0, ddim)
  order[mask] <- z$order
  lev <- matrix(0,2,nvox)
  lev[2,mask] <- log(z$lambda)
  lev[1,mask] <- log((1+z$alpha)*z$lambda)
  dim(lev) <- c(2,ddim)
  mix <- matrix(0,maxcomp,nvox)
  mix[,mask] <- z$mix
  dim(mix) <- c(maxcomp, ddim)
} else {
  cat("Starting parameter estimation and model selection (C-code) on",mc.cores," cores",format(Sys.time()),"\n")
  x <- matrix(0,ngrad0+3+maxcomp,sum(mask))
  dim(siq) <- c(ngrad0,nvox)
  x[1:ngrad0,] <- siq[,mask]
  x[ngrad0+1,] <- sigma2[mask]
  dim(siind) <- c(2+maxcomp,nvox)
  x[ngrad0+2:(1+maxcomp),] <- siind[,mask] 
  res <- matrix(0,4+3*maxcomp,nvox)
  res[,mask] <- switch(model+1,
                  plmatrix(x,pmixtns0,ngrad0=ngrad0,maxcomp=maxcomp,maxit=maxit,
                      grad=grad,bv=bvalue,lambda=lambda,alpha=alpha,factr=factr,
                      penIC=penIC,vert=vert,mc.cores=mc.cores),# model=0
                  plmatrix(x,pmixtns1,ngrad0=ngrad0,maxcomp=maxcomp,maxit=maxit,
                      grad=grad,bv=bvalue,lambda=lambdahat,alpha=alpha,factr=factr,
                      penIC=penIC,vert=vert,mc.cores=mc.cores),# model=1
                  plmatrix(x,pmixtns1,ngrad0=ngrad0,maxcomp=maxcomp,maxit=maxit,
                      grad=grad,bv=bvalue,lambda=lambdahat,alpha=alphahat,factr=factr,
                      penIC=penIC,vert=vert,mc.cores=mc.cores))# model=2
  cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
  rm(x)
  gc()
  sigma2 <-  array(res[2,],ddim)
  orient <- array(res[maxcomp+4+1:(2*maxcomp),], c(2, maxcomp, ddim))
  order <- array(as.integer(res[1,]), ddim)
  lev <- array(res[3:4,], c(2,ddim))
  mix <- array(res[4+(1:maxcomp),], c(maxcomp, ddim))
  }
  method <- switch(model,"mixtensoriso0","mixtensoriso1","mixtensoriso2")
  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
                call   = args,
                ev     = lev-log(maxbv),# reverse scaling of eigenvalues
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
)
#
#   create initial estimates for directions
#
getsiind2 <- function(si,mask,sigma2,grad,bv,vico,alpha,lambda,maxcomp=3,maxc=.866,nguess=100){
# assumes dim(grad) == c(ngrad,3)
# assumes dim(si) == c(ngrad,n1,n2,n3)
# SO removed
ngrad <- dim(grad)[1]
nvico <- dim(vico)[1]
nsi <- dim(si)[1]
dgrad <- matrix(abs(grad%*%t(vico)),ngrad,nvico)
dgrad <- dgrad/max(dgrad)
dgradi <- matrix(abs(vico%*%t(vico)),nvico,nvico)
dgradi <- dgradi/max(dgradi)
isample <- selisample(nvico,maxcomp,nguess,dgradi,maxc)
#
#  eliminate configurations with close directions 
#
# this provides configurations of initial estimates with minimum angle between 
# directions > acos(maxc)
nvoxel <- prod(dim(si)[-1])
cat("using ",nguess,"guesses for initial estimates\n")
siind <- matrix(as.integer(0),maxcomp+1,nvoxel)
krit <- numeric(nvoxel)
nguess <- length(isample)/maxcomp
z <- .Fortran("getsii",
         as.double(si),
         as.double(sigma2),
         as.integer(nsi),
         as.integer(nvoxel),
         as.integer(maxcomp),
         as.double(dgrad),
         as.double(bv),
         as.integer(nvico),
         as.double(alpha),
         as.double(lambda),
         double(ngrad*nvico),
         as.integer(isample),
         as.integer(nguess),
         double(nsi),
         double(nsi),#z0
         double(nsi*(maxcomp+1)),
         siind=integer((maxcomp+1)*nvoxel),
         krit=double(nvoxel),
         as.logical(mask),
         as.integer(maxcomp+1),
         PACKAGE="dti")[c("siind","krit")]
dim(z$siind) <- c(maxcomp+1,nvoxel)
siind <- z$siind
krit <- z$krit
# now voxel where first tensor direction seems important
failed <- (krit^2/ngrad) > (sigma2-1e-10)
#if(any(failed[mask])){
#print((krit[mask])[failed[mask]])
#print(((1:prod(dim(si)[1:3]))[mask])[failed[mask]])
#print(sum(failed[mask]))
#}
list(siind=array(siind,c(maxcomp+1,dim(si)[-1])),
     krit=array(krit,dim(si)[-4]))
}
mixtens <-
function(object, maxcomp=3, fa=NULL, lambda=NULL, factr=1e7, maxit=5000,ngc=1000, nguess=100*maxcomp^2,msc="BIC",thinit=NULL, 
    mc.cores = setCores(,reprt=FALSE)){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
#  uses L-BFGS-B  
#  models with isotropic compartment
#     if(!is.null(FA)&!is.null(lambda)) use max. EV / FA as provided
#     if(!is.null(FA)&is.null(lambda)) use FA as provided
#     if(is.null(FA)) estimate both EV's
#
  set.seed(1)
#
#  pre-processing
#
  theta <- .5
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  nvox <- prod(ddim)
  s0ind <- object@s0ind
  bvalue <- object@bvalue[-s0ind]
  maxbv <- max(bvalue)
  bvalue <- bvalue/maxbv
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(5*(1+3*maxcomp)>ngrad0){
#     maxcomp <- max(1,trunc((ngrad0-5)/15))
     cat("Maximal number of components reduced to", maxcomp,"due to insufficient
          number of gradient directions\n")
  }
#
#   which model should be used
#
if(is.null(fa)){
   model <- 2  
} else {
   fa2 <- fa^2
   alpha <- (fa2+sqrt(fa2*(3-2*fa2)))/(1-fa2)
   if(is.null(lambda)){
      model <- 1
   } else {
      model <- 0
      lambda <- lambda*maxbv/(1+alpha)
   }
}
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
  lambdahat <- ev[3,,,] 
# use third ev instead of (ev[2,,,]+ev[3,,,])/2 to avoid effects from mixtures
  alphahat <- if(model==2) (ev[1,,,]-lambdahat)/lambdahat else alpha
  lambdahat <- median(lambdahat[!is.na(alphahat)&fa>.3])
  if(model==2) alphahat <- median(alphahat[!is.na(alphahat)&fa>.3])
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
  npar <- model+3*(0:maxcomp)
#
#   compute penalty for model selection, default BIC
#
  penIC <- switch(msc,"AIC"=2*npar/ngrad0,"BIC"=log(ngrad0)*npar/ngrad0,
                  "AICC"=(1+npar/ngrad0)/(1-(npar+2)/ngrad0),
                  "None"=log(ngrad0)-log(ngrad0-npar),
                  log(ngrad0)*npar/ngrad0)
  cat("End generating auxiliary objects",format(Sys.time()),"\n")
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
#
#   determine initial estimates for orientations 
#
  cat("Start search for initial directions at",format(Sys.time()),"\n")
  data("polyeders")
  polyeder <- icosa3
  vert <- polyeder$vertices
# remove redundant directions
  vind <- rep(TRUE,dim(vert)[2])
  vind[vert[1,]<0] <- FALSE
  vind[vert[1,]==0 & vert[2,] <0] <- FALSE
  vind[vert[1,]==0 & vert[2,] == 0 &vert[3,]<0] <- FALSE
  vert <- vert[,vind]
#
#  compute initial estimates (EV from grid and orientations from icosa3$vertices)
#
  siind <-   getsiind2(siq,mask,sigma2,grad,bvalue,t(vert),alphahat,lambdahat,
                       maxcomp,maxc=maxc,nguess=nguess)
  krit <- siind$krit # sqrt(sum of squared residuals) for initial estimates
  siind <- siind$siind # components 1: model order 2: 
                       # grid index for EV 2+(1:m) index of orientations
  cat("Model orders for initial estimates")
  print(table(siind[1,,,]))
  cat("End search for initial values at",format(Sys.time()),"\n")
  siind <- siind[-1,,,]
# only need indices of initial directions
  orient <- array(0,c(2,maxcomp,ddim))
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  prt0 <- Sys.time()
#
#     C-Code
#
if(mc.cores<=1){
  cat("Starting parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
  dim(siq) <- c(ngrad0,nvox)
  dim(siind) <- c(maxcomp,nvox)
  nvoxm <- sum(mask)
  z <- switch(model+1,.C("mixtrl0", 
          as.integer(nvoxm),#n1
          as.integer(siind[,mask]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bvalue),#bv_in
          lambda  = as.double(lambda),#lambda_in
          alpha   = as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(sigma2[mask]),#sigma2
          as.double(vert),#vert
          as.double(siq[,mask]),#siq_in
          sigma2  = double(nvoxm),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvoxm),#order_ret selected order of mixture
          mix     = double(maxcomp*nvoxm),#mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")],
          .C("mixtrl1", 
          as.integer(nvoxm),#n1
          as.integer(siind[,mask]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bvalue),#bv_in
          as.double(lambdahat),#lambda_in
          alpha = as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(sigma2[mask]),#sigma2
          as.double(vert),#vert
          as.double(siq[,mask]),#siq_in
          sigma2  = double(nvoxm),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvoxm),#order_ret selected order of mixture
          lambda  = double(nvoxm),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvoxm),#mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")],
          .C("mixtrl2", 
          as.integer(nvoxm),#n1
          as.integer(siind[,mask]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bvalue),#bv_in
          as.double(lambdahat),#lambda_in
          as.double(alphahat),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(sigma2[mask]),#sigma2
          as.double(vert),#vert
          as.double(siq[,mask]),#siq_in
          sigma2  = double(nvoxm),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvoxm),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvoxm),#order_ret selected order of mixture
          alpha   = double(nvoxm),#alpha_ret alpha=(lambda_1-lambda_2)/lambda_2 
          lambda  = double(nvoxm),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvoxm),#mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")])
  cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
  sigma2 <-  array(0,ddim)
  sigma2[mask] <- z$sigma2
  orient <- matrix(0,2*maxcomp,nvox)
  orient[,mask] <- z$orient
  dim(orient) <- c(2, maxcomp, ddim)
  order <- array(0, ddim)
  order[mask] <- z$order
  lev <- matrix(0,2,nvox)
  lev[2,mask] <- z$lambda
  lev[1,mask] <- z$alpha*z$lambda
  dim(lev) <- c(2,ddim)
  mix <- matrix(0,maxcomp,nvox)
  mix[,mask] <- z$mix
  dim(mix) <- c(maxcomp, ddim)
} else {
  cat("Starting parameter estimation and model selection (C-code) on",mc.cores," cores",format(Sys.time()),"\n")
  x <- matrix(0,ngrad0+3+maxcomp,sum(mask))
  dim(siq) <- c(ngrad0,nvox)
  x[1:ngrad0,] <- siq[,mask]
  x[ngrad0+1,] <- sigma2[mask]
  dim(siind) <- c(maxcomp,nvox)
  x[ngrad0+2:(1+maxcomp),] <- siind[,mask] 
  res <- matrix(0,4+3*maxcomp,nvox)
  res[,mask] <- switch(model+1,
                  plmatrix(x,pmixtns0,ngrad0=ngrad0,maxcomp=maxcomp,maxit=maxit,
                      grad=grad,bv=bvalue,lambda=lambda,alpha=alpha,factr=factr,
                      penIC=penIC,vert=vert,mc.cores=mc.cores),# model=0
                  plmatrix(x,pmixtns1,ngrad0=ngrad0,maxcomp=maxcomp,maxit=maxit,
                      grad=grad,bv=bvalue,lambda=lambdahat,alpha=alpha,factr=factr,
                      penIC=penIC,vert=vert,mc.cores=mc.cores),# model=1
                  plmatrix(x,pmixtns1,ngrad0=ngrad0,maxcomp=maxcomp,maxit=maxit,
                      grad=grad,bv=bvalue,lambda=lambdahat,alpha=alphahat,factr=factr,
                      penIC=penIC,vert=vert,mc.cores=mc.cores))# model=2
  cat("End parameter estimation and model selection (C-code)",format(Sys.time()),"\n")
  rm(x)
  gc()
  sigma2 <-  array(res[2,],ddim)
  orient <- array(res[maxcomp+4+1:(2*maxcomp),], c(2, maxcomp, ddim))
  order <- array(as.integer(res[1,]), ddim)
  lev <- array(res[3:4,], c(2,ddim))
  mix <- array(res[4+(1:maxcomp),], c(maxcomp, ddim))
  }
  method <- switch(model,"mixtensoriso0","mixtensoriso1","mixtensoriso2")
  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
                call   = args,
                ev     = lev/maxbv,# reverse scaling of eigenvalues
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

