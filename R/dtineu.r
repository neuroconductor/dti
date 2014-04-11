################################################################
#                                                              #
# Section for DTI functions                                    #
#                                                              #
################################################################

dtiTensor <- function(object,  ...) cat("No DTI tensor calculation defined for this class:",class(object),"\n")

setGeneric("dtiTensor", function(object,  ...) standardGeneric("dtiTensor"))

setMethod( "dtiTensor", "dtiData",
           function( object, 
                     method = c( "nonlinear", "linear"),
                     mc.cores = setCores( , reprt = FALSE)) {

## check method! available are: 
##   "linear" - use linearized model (log-transformed)
##   "nonlinear" - use nonlinear model with parametrization according to Koay et.al. (2006)
             method <- match.arg( method)
                          
   if(is.null(mc.cores)) mc.cores <- 1
   mc.cores <- min(mc.cores,detectCores())
   args <- sys.call(-1)
   args <- c(object@call,args)
   ngrad <- object@ngrad
   grad <- object@gradient
   ddim <- object@ddim
   ntotal <- prod(ddim)
   s0ind <- object@s0ind
   ns0 <- length(s0ind)
   sdcoef <- object@sdcoef
   ngrad0 <- ngrad-ns0
   #if(mc.cores>1&method=="nonlinear") require(parallel)
   if(all(sdcoef[1:4]==0)) {
      cat("No parameters for model of error standard deviation found\n
          Using constant weights \n You may prefer to run sdpar before calling dtiTensor\n")
#      sdcoef <- sdpar(object,interactive=FALSE)@sdcoef
       sdcoef <- c(1,0,1,1)
   }
   z <- sioutlier1(object@si,s0ind,object@level,mc.cores=mc.cores)
#
#  this does not scale well with openMP
#
cat("sioutlier completed\n")
   mask <- z$mask
   nvox <- sum(mask)
   ttt <- array(0,c(ngrad0,nvox))
   ttt <- -log1p(sweep(z$si[-s0ind,],2,as.vector(z$s0),"/")-1)
#  suggestion by B. Ripley
#   idsi <- 1:length(dim(si))
#   ttt <- -log(sweep(si,idsi[-1],s0,"/"))
   ttt[is.na(ttt)] <- 0
   ttt[(ttt == Inf)] <- 0
   ttt[(ttt == -Inf)] <- 0
   cat("Data transformation completed ",format(Sys.time()),"\n")
               
   D <- matrix(0,6,ntotal)
   btbsvd <- svd(object@btb[,-s0ind])
   solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
   D[,mask] <- solvebtb%*% ttt
   cat("Diffusion tensors (linearized model) generated ",format(Sys.time()),"\n")
   rm(ttt)
   D[c(1,4,6),!mask] <- 1e-6
   D[c(2,3,5),!mask] <- 0
   D <- dti3Dreg(D,mc.cores=mc.cores)
   dim(D) <- c(6,ntotal)
   th0 <- array(0,ntotal)
   th0[mask] <- z$s0
   index <- z$index
   if(method== "nonlinear"){
#  method == "nonlinear" uses linearized model for initial estimates
       ngrad0 <- ngrad
       df <- sum(table(object@replind)-1)
       param <- matrix(0,7,nvox)
       ms0 <- mean(z$s0)
       mbv <- mean(object@bvalue[-object@s0ind])
##  use ms0 and mbv to rescale parameters such that they are of comparable magnitude
       param[1,] <- z$s0/ms0
       param[-1,] <- D2Rall(D[,mask]*mbv)
## use reparametrization D = R^T R
       sdcoef[-2] <- sdcoef[-2]/ms0# effect of rescaling of signal
       cat("start nonlinear regression",format(Sys.time()),"\n")
       if(mc.cores==1){
          param <- matrix(.C("dtens",
                         as.integer(nvox),
                         param=as.double(param),
                         as.double(z$si/ms0),
                         as.integer(ngrad),
                         as.double(object@btb/mbv),
                         as.double(sdcoef),
                         as.double(rep(0,ngrad)),#si
                         as.double(rep(1,ngrad)),#var                         
                         as.integer(1000),#maxit
                         as.double(1e-7),#reltol
                         PACKAGE="dti",DUP=TRUE)$param,7,nvox)
       } else {
          x <- matrix(0,ngrad+7,nvox)
          x[1:7,] <- param
          x[-(1:7),] <- z$si/ms0
          param <- plmatrix(x,ptensnl,ngrad=ngrad,btb=object@btb/mbv,
                            sdcoef=sdcoef,maxit=1000,reltol=1e-7)
       }
       th0[mask] <- param[1,]*ms0
       D[,mask] <- R2Dall(param[-1,])/mbv
       cat("successfully completed nonlinear regression ",format(Sys.time()),"\n")
   }
   res <- matrix(0,ngrad,ntotal)
   z <- .Fortran("tensres",
                 as.double(th0[mask]),
                 as.double(D[,mask]),
                 as.double(z$si),
                 as.integer(nvox),
                 as.integer(ngrad),
                 as.double(object@btb),
                 res = double(ngrad*nvox),
                 rss = double(nvox),
                 PACKAGE="dti",DUP=TRUE)[c("res","rss")]
    res <- matrix(0,ngrad,ntotal)
    res[,mask] <- z$res
    rss <- numeric(ntotal)
    rss[mask] <- z$rss/(ngrad0-6)
    rm(z)
    dim(th0) <- ddim
    dim(D) <- c(6,ddim)
    dim(res) <- c(ngrad,ddim)
    dim(rss) <- ddim
    dim(mask) <-  ddim
    gc()
#
#   get spatial correlation
#
    if(any(is.na(res))){
       dim(res) <- c(ngrad,ntotal)
       indr <- (1:ntotal)[apply(is.na(res),2,any)]
       cat("NA's in res in voxel",indr,"\n")
       res[,indr] <- 0
    }
    if(any(is.na(D))|any(abs(D)>1e10)){
       dim(D) <- c(6,total)
       indD <- (1:total)[apply(is.na(D),2,any)]
       cat("NA's in D in ", length(indD),"voxel:",indD,"\n")
       D[,indD] <- c(1,0,0,1,0,1)
       mask[indD] <- FALSE
       indD <- (1:ntotal)[apply(abs(D)>1e10,2,any)]
       cat("Inf's in D in", length(indD)," voxel:",indD,"\n")
       D[,indD] <- c(1,0,0,1,0,1)
       mask[indD] <- FALSE
    }
    scorr <- mcorr(res,mask,ddim,ngrad0,lags=c(5,5,3),mc.cores=mc.cores)
    ev <- dti3Dev(D,mask,mc.cores=mc.cores)
    dim(ev) <- c(3,ddim)   
    dim(D) <- c(6,ddim)   
    scale <- quantile(ev[3,,,][mask],.95,na.rm=TRUE)
    cat("estimated scale information",format(Sys.time()),"\n")  
    invisible(new("dtiTensor",
                  call  = args,
                  D     = D,
                  th0   = th0,
                  sigma = rss,
                  scorr = scorr$scorr, 
                  bw = scorr$bw, 
                  mask = mask,
                  hmax = 1,
                  gradient = object@gradient,
                  bvalue = object@bvalue,
                  btb   = object@btb,
                  ngrad = object@ngrad, # = dim(btb)[2]
                  s0ind = object@s0ind,
                  replind = object@replind,
                  ddim  = object@ddim,
                  ddim0 = object@ddim0,
                  xind  = object@xind,
                  yind  = object@yind,
                  zind  = object@zind,
                  voxelext = object@voxelext,
                  level = object@level,
                  orientation = object@orientation,
                  rotation = object@rotation,
                  source = object@source,
                  outlier = index,
                  scale = scale,
                  method = method)
             )
})

D2Rall <- function(D){
#
#  used for initial values
# in case of negative eigenvalues this uses c(1,0,0,1,0,1) for initial rho
#
      nvox <- dim(D)[2]
      matrix(.Fortran("D2Rall",
               as.double(D),
               rho=double(6*nvox),
               as.integer(nvox),
               DUP=FALSE,
               PACKAGE="dti")$rho,6,nvox)
}
R2Dall <- function(R){
#
#  used for initial values
# in case of negative eigenvalues this uses c(1,0,0,1,0,1) for initial rho
#
      nvox <- dim(R)[2]
      matrix(.Fortran("R2Dall",
               as.double(R),
               D=double(6*nvox),
               as.integer(nvox),
               DUP=FALSE,
               PACKAGE="dti")$D,6,nvox)
}
#############
#############

dtiIndices <- function(object, ...) cat("No DTI indices calculation defined for this class:",class(object),"\n")

setGeneric("dtiIndices", function(object, ...) standardGeneric("dtiIndices"))

setMethod("dtiIndices","dtiTensor",
function(object, mc.cores=setCores(,reprt=FALSE)) {
  args <- sys.call(-1)
  args <- c(object@call,args)
  ddim <- object@ddim
  n <- prod(ddim)
  z <- dtiind3D(object@D,object@mask,mc.cores=mc.cores)
  invisible(new("dtiIndices",
                call = args,
                fa = array(z$fa,ddim),
                ga = array(z$ga,ddim),
                md = array(z$md,ddim),
                andir = array(z$andir,c(3,ddim)),
                bary = array(z$bary,c(3,ddim)),
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = ddim,
                ddim0 = object@ddim0,
                voxelext = object@voxelext,
                orientation = object@orientation,
                rotation = object@rotation,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                method = object@method,
                level = object@level,
                source= object@source)
            )
})

dtiTensorChi <- function(object, sigma=NULL, L=1,
                         mc.cores = setCores( , reprt = FALSE)){
##
##  Estimate parameters in Diffusion Tensor model with
##  Gauss-approximation for noncentral chi
##
   tensobj <- dtiTensor(object,method = "nonlinear")
   if(is.null(sigma)||sigma<=0){
      cat("Please specify a valid sigma\n Returning biased estimated tensor object\n")
      return(tensobj)
   }
   ddim <- tensobj@ddim
   bv <- object@bvalue
   D <- tensobj@D
   dim(D) <- c(6,prod(ddim))
   mbv <- mean(bv[-object@s0ind])
   btb <- object@btb/mbv
##  use ms0 and mbv to rescale parameters such that they are of comparable magnitude
   z <- sioutlier1(object@si,object@s0ind,object@level,mc.cores=mc.cores)
   mask <- z$mask
   nvox <- sum(mask)
   
   param <- matrix(0,7,nvox)
   if(length(sigma)==1) sigma <- array(sigma,ddim)
   param[1,] <- tensobj@th0[mask]/sigma[mask]
   param[-1,] <- D2Rall(D[,mask]*mbv)
   
   CL <- sqrt(pi/2)*gamma(L+1/2)/gamma(L)/gamma(3/2)
   if(mc.cores==1){
      si <- t(array(z$si,c(object@ngrad,nvox)))/sigma[mask]
      for(i in 1:length(mask)){
         param[,i] <- optim(param[,i],tchi,si=si[,i],btb=btb,L=L,CL=CL,method="BFGS",
                            control=list(reltol=1e-5,maxit=50))$par
         if(i%/%1000*1000==i) cat(i,"voxel processed. Time:",format(Sys.time()),"\n")
   }
   } else {
      x <- matrix(0,object@ngrad+7,nvox)
      x[1:7,] <- param
      x[-(1:7),] <- t(t(array(z$si,c(object@ngrad,nvox)))/sigma[mask])
      param <- plmatrix(x,ptenschi,fn=tchi,btb=btb,L=L,CL=CL)
      cat(nvox,"voxel processed. Time:",format(Sys.time()),"\n")
   }
   D[,mask] <- R2Dall(param[-1,])/mbv
   dim(D) <- c(6,ddim)
   tensobj@D <- D
   tensobj@th0[mask] <- param[1,]*sigma[mask]
   tensobj
}
##
## Parallel version
##
ptenschi <- function(x,fn,btb,L,CL){
nvox <- dim(x)[2]
param <- matrix(0,7,nvox)
for(i in 1:nvox){
   param[,i] <-
   optim(x[1:7,i],fn,si=x[-(1:7),i],btb=btb,L=L,CL=CL,method="BFGS",
                      control=list(reltol=1e-5,maxit=50))$par
   }
param
}


tchi <- function(param,si,btb,L,CL){
##
##  Risk function for Diffusion Tensor model with
##  Gauss-approximation for noncentral chi
##
   ng <- dim(btb)[2]
   D <- .Fortran("rho2D0",
                 as.double(param[-1]),
                 D=double(6),
                 DUPL=FALSE,
                 PACKAGE="dti")$D
   gDg <- D%*%btb ## b_i*g_i^TD g_i (i=1,ngrad)
   gvalue <- param[1]*exp(-gDg)
   mgvh <- -gvalue*gvalue/2
   muL <- CL*.C("hyperg_1F1_e",
              as.double(rep(-.5,ng)),
              as.double(rep(L,ng)),
              as.double(mgvh),
              as.integer(ng),
              val=as.double(mgvh),
              err=as.double(mgvh),
              status=as.integer(0*mgvh),
              PACKAGE="gsl")$val
   vL <- 2*L+gvalue^2-muL^2
   sum((si-muL)^2/vL)
}

