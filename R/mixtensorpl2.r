dwiMixtensor2 <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor2", function(object,  ...) standardGeneric("dwiMixtensor2"))

setMethod("dwiMixtensor2","dtiData",function(object, maxcomp=3, method="mixtensor", reltol=1e-6, maxit=5000,ngc=1000, optmethod="BFGS", nguess=100*maxcomp^2,msc="BIC",pen=NULL){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
#  choices for optmethod:
#  BFGS  -  BFGS with analytic gradients and penalization
#  CG - Conjugate gradients with analytic gradients and penalization
#  L-BFGS-B  -  constrained BFGS with analytic gradients 
#  Nelder-Mead - using LawsonHanson-nnls code
#
#  Defaults: 
#     BFGS for tensor mixture models without isotropic compartment
#     L-BFGS-B for tensor mixture models with isotropic compartment
#
  set.seed(1)
  if(is.null(pen)) pen <- 100
  if(method=="mixtensoriso") optmethod <- "L-BFGS-B"
  theta <- .5
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(5*(1+3*maxcomp)>ngrad0){
#     maxcomp <- max(1,trunc((ngrad0-5)/15))
     cat("Maximal number of components reduced to", maxcomp,"due to insufficient
          number of gradient directions\n")
  }


#
#  First tensor estimates to generate eigenvalues and -vectors
#
  prta <- Sys.time()
  cat("Start tensor estimation at",format(prta),"\n")
  tensorobj <- dtiTensor(object)
  cat("Start evaluation of eigenstructure at",format(Sys.time()),"\n")
  z <- .Fortran("dtieigen",
                as.double(tensorobj@D),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.logical(tensorobj@mask),
                fa=double(prod(ddim)),
                ev=double(3*prod(ddim)),
                andir=double(6*prod(ddim)),
                DUP=FALSE,
                PACKAGE="dti")[c("fa","ev","andir")]
  rm(tensorobj)
  gc()
  fa <- array(z$fa,ddim)
  ev <- array(z$ev,c(3,ddim))
  andir <- array(z$andir,c(3,2,ddim))
  rm(z)
  gc()
  nth <- 11
  th <- ev[1,,,] - (ev[2,,,]+ev[3,,,])/2
  falevel <- min(quantile(fa[fa>0],.75),.3)
  rth <- quantile(th[fa>=falevel],c(.1,.99))
  if(diff(rth)>0){
     indth <- trunc((nth-1)*(th-rth[1])/diff(rth)+1)
     th <- seq(rth[1],rth[2],length=nth)
  } else {
     th <- rep(max(th),nth)
     indth <- rep(1,length(th))
  }
cat("using th:::",th,"\n")
  cat("Start search outlier detection at",format(Sys.time()),"\n")

#
#  replace physically meaningless S_i by mena S_0 values
#
  z <- .Fortran("outlier",#misc.f 
                as.double(object@si),
                as.integer(prod(ddim)),
                as.integer(ngrad),
                as.logical((1:ngrad)%in%s0ind),
                as.integer(ns0),
                si=integer(prod(ddim)*ngrad),
                index=integer(prod(ddim)),
                lindex=integer(1),
                DUP=FALSE,
                PACKAGE="dti")[c("si","index","lindex")]
  cat("End search outlier detection at",format(Sys.time()),"\n")
  si <- array(as.integer(z$si),c(ddim,ngrad))
  index <- if(z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  gc()
  cat("Start generating auxiliary objects",format(Sys.time()),"\n")

#
#  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
#
  z <- .Fortran("sweeps0",# mixtens.f
                as.integer(si[,,,-s0ind,drop=FALSE]),
                as.integer(si[,,,s0ind,drop=FALSE]),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.integer(ns0),
                as.integer(ngrad0),
                as.integer(object@level),
                siq=double(prod(ddim[1:3])*ngrad0),
                s0=double(prod(ddim[1:3])),
                vsi=double(prod(ddim[1:3])),
                mask=logical(prod(ddim[1:3])),
                DUPL=FALSE,
                PACKAGE="dti")[c("siq","s0","vsi","mask")]
  rm(si)
  s0 <- array(z$s0,ddim[1:3])
  siq <- array(z$siq,c(ddim[1:3],ngrad0))
  sigma2 <- array(z$vsi,ddim[1:3])
  mask <- array(z$mask,ddim[1:3])
  rm(z)
  gc()
  npar <- if(method=="mixtensor") 1+3*(0:maxcomp) else c(1,2+3*(1:maxcomp))

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
  siind <- if(method=="mixtensor")  
    getsiind3(siq,mask,sigma2,grad,t(vert),th,indth,ev,fa,andir,maxcomp,maxc=maxc,nguess=nguess) else 
    getsiind3iso(siq,mask,sigma2,grad,t(vert),th,indth,ev,fa,andir,maxcomp,maxc=maxc,nguess=nguess)

  krit <- siind$krit # sqrt(sum of squared residuals) for initial estimates
  siind <- siind$siind # components 1: model order 2: 
                       # grid index for EV 2+(1:m) index of orientations
  cat("Model orders for initial estimates")
  print(table(siind[1,,,]))
  cat("End search for initial values at",format(Sys.time()),"\n")
  order <- array(0,ddim)

#  logarithmic eigen values
  mix <- array(0,c(maxcomp,ddim))
  orient <- array(0,c(2,maxcomp,ddim))
  lev <- array(0,c(2,ddim))
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  prt0 <- Sys.time()

  if(method=="mixtensor") meth = 1 else meth = 2
  optmeth <- switch(optmethod,   "BFGS" = 1,
			    "CG"   = 2,
			    "Nelder-Mead" = 3,
			    "L-BFGS-B" = 4)


  cat("initial parameters: ", n1,n2,n3,ngrad,ngrad0, "\n")
  print("starting C code...")
  print("method")
  print(optmeth)


  z <- .C("mixture2", 
          as.integer(meth),
          as.integer(optmeth), 
          as.integer(n1), 
          as.integer(n2), 
          as.integer(n3),
          as.integer(mask), 
          as.integer(siind), 
          as.integer(ngrad), 
          as.integer(ngrad0),
          as.integer(maxcomp),
          as.integer(maxit),
          as.double(pen),
          as.double(t(grad)),
          as.double(reltol),
          as.double(th),
          as.double(penIC),
          as.double(sigma2),
          as.double(vert),
          as.double(orient),
          as.double(siq),
          sigma2  = double(prod(ddim)),           # parameter variance ???
          orient  = double(2*maxcomp*prod(ddim)), # phi/theta for all mixture tensors
          order   = integer(prod(ddim)),   # selected order of mixture
          lev     = double(2*prod(ddim)),         # logarithmic eigenvalues
          mix     = double(maxcomp*prod(ddim)),   # mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","vert","orient","order","lev","mix")]

sigma2 <-  array(z$sigma2,ddim[1:3])
orient <- array(z$orient, c(2, maxcomp, ddim[1:3]))
order <- array(z$order, ddim[1:3])
lev <- array(z$lev, c(2,ddim[1:3]))
mix <- array(z$mix, c(maxcomp, ddim[1:3]))
method <- "mixtensor"


  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
                call   = args,
                ev     = lev,
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

 