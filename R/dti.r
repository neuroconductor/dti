  # solvebtb mit einfuegen??

setMethod("show", "dti",
function(object){
    cat("DTI object\n")
    cat("  Dimension     :", paste(object@ddim, collapse="x"), "\n")
    cat("  Filename      :", object@source, "\n")
    cat("\n")
})

setMethod("plot", "dtiTensor", function(x, y, ...) cat("Not yet implemented yet for class dtiTensor\n"))
setMethod("plot", "dtiData", function(x, y, ...) cat("Not yet implemented yet for class dtiData\n"))
setMethod("plot", "dti", function(x, y, ...) cat("No implementation for class dti\n"))

setMethod("plot", "dtiIndices", 
function(x, y, slice=1, method=1, quant=0, minanindex=NULL, show=TRUE, fa.exp = 1, ...) {
  if(is.null(x$fa)) cat("No anisotropy index yet")
  adimpro <- require(adimpro)
  anindex <- x$fa[,,slice]
  dimg <- x@ddim[1:2]
  if ((method==1) || (method==2)) {
    andirection <- x$eigenv[,,slice,,]
    anindex[anindex>1]<-0
    anindex[anindex<0]<-0
    if ((fa.exp != 1) & (fa.exp>0)) anindex <- anindex^fa.exp
    dim(andirection)<-c(prod(dimg),3,3)
    if(is.null(minanindex)) minanindex <- quantile(anindex,quant,na.rm=TRUE)
    if (diff(range(anindex,na.rm=TRUE)) == 0) minanindex <- 0
    if(method==1) {
      andirection[,1,3] <- andirection[,1,3]^2
      andirection[,2,3] <- andirection[,2,3]^2
      andirection[,3,3] <- andirection[,3,3]^2
    } else {
      ind<-andirection[,1,3]<0
      andirection[ind,,] <- - andirection[ind,,]
      andirection[,2,3] <- (1+andirection[,2,3])/2
      andirection[,3,3] <- (1+andirection[,3,3])/2
    }
    andirection <- andirection[,,3]
    andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
    dim(andirection)<-c(dimg,3)
    if(adimpro) {
      andirection[is.na(andirection)] <- 0
      andirection <- make.image(andirection)
      if(show) show.image(andirection,...)
    } else if(show) {
      dim(anindex) <- dimg
      image(anindex,...)
    }
    invisible(andirection)
  } else if (method==3) {
    bary <- x$bary[,,slice,]
    if(adimpro) {
      bary[is.na(bary)] <- 0
      bary <- make.image(bary)
      if(show) show.image(bary,...)
    } else if(show) {
      image(bary[,,1],...)
    }
    invisible(bary)
  } else {
    andirection <- x$eigenv[,,slice,,]
    anindex[anindex>1]<-0
    anindex[anindex<0]<-0
    if ((fa.exp != 1) & (fa.exp>0)) anindex <- anindex^fa.exp
    dim(andirection)<-c(prod(dimg),3,3)
    if(is.null(minanindex)) minanindex <- quantile(anindex,quant,na.rm=TRUE)
    if (diff(range(anindex,na.rm=TRUE)) == 0) minanindex <- 0
    andirection[,1,2] <- andirection[,1,2]^2
    andirection[,2,2] <- andirection[,2,2]^2
    andirection[,3,2] <- andirection[,3,2]^2
    andirection <- andirection[,,2]
    andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
    dim(andirection)<-c(dimg,3)
    if(adimpro) {
      andirection[is.na(andirection)] <- 0
      andirection <- make.image(andirection)
      if(show) show.image(andirection,...)
    } else if(show) {
      dim(anindex) <- dimg
      image(anindex,...)
    }
    invisible(andirection)
  }
})


dtiData <- function(gradient,imagefile,ddim,xind=NULL,yind=NULL,zind=NULL,level=0,voxelext=c(1,1,1)) {
  if (dim(gradient)[2]==3) gradient <- t(gradient)
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]

  if (!(file.exists(imagefile))) stop("Image file does not exist")
  zz <- file(imagefile,"rb")
  s0 <- readBin(zz,"integer",prod(ddim),2,FALSE)
  si <- readBin(zz,"integer",prod(ddim)*ngrad,2,FALSE)
  close(zz)
  cat("Data successfully read \n")

  if (is.null(xind)) xind <- 1:ddim[1]
  if (is.null(yind)) yind <- 1:ddim[2]
  if (is.null(zind)) zind <- 1:ddim[3]
  dim(s0) <- ddim
  s0 <- s0[xind,yind,zind] # really needed?
  dim(si) <- c(ddim,ngrad)
  si <- si[xind,yind,zind,] # really needed?
  ddim0 <- as.integer(ddim)
  ddim <- dim(s0)

  btb <- create.designmatrix.dti(gradient)

  invisible(new("dtiData",
                list(s0 = s0, si = si),
                btb    = btb,
                ngrad  = ngrad, # = dim(btb)[2]
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = xind,
                yind   = yind,
                zind   = zind,
                level  = level,
                voxelext = voxelext,
                source = imagefile)
            )
}

# has to be re-implemented!!!!!!!!!!!!!!!!!!!!!!!!!!1
createdata.dti <- function(file,dtensor,btb,s0,sigma,level=250){
  ngrad <- dim(btb)[2]
  ddim <- dim(s0)
  dtensor[1,,,][s0<level] <- 1e-5
  dtensor[4,,,][s0<level] <- 1e-5
  dtensor[6,,,][s0<level] <- 1e-5
  dtensor[2,,,][s0<level] <- 0
  dtensor[3,,,][s0<level] <- 0
  dtensor[5,,,][s0<level] <- 0
  dim(dtensor)<-c(6,prod(ddim))
  dtensor <- t(dtensor)
  si <- exp(-dtensor%*%btb)*as.vector(s0)
  rsi <- pmax(0,rnorm(si,si,pmin(s0/2.5,sigma)))
  rs0 <- pmax(0,rnorm(s0,s0,pmin(s0/2.5,sigma)))
  zz <- file(file,"wb")
  writeBin(as.integer(rs0),zz,2)
  writeBin(as.integer(rsi),zz,2)
  close(zz)
  dim(s0)<-ddim
  dim(si)<-c(ddim,ngrad)
  dtensor <- t(dtensor)
  dim(dtensor)<-c(6,ddim)
  list(s0=s0,si=si,dtensor=dtensor,sigma=sigma,level=level,btb=btb)
}
# END implementation needed!

# really setAs() or setMethod?
setAs("dtiData","dtiTensor",function(from,to) {
  ngrad <- from@ngrad
  ddim <- from@ddim

  s0 <- from$s0
  si <- from$si
  dim(s0) <- dim(si) <- NULL
  ttt <- -log(si/s0)
  ttt[is.na(ttt)] <- 0
  ttt[(ttt == Inf)] <- 0
  ttt[(ttt == -Inf)] <- 0
  dim(ttt) <- c(prod(ddim),ngrad)
  ttt <- t(ttt)
  cat("Data transformation completed \n")

  btbsvd <- svd(from@btb)
  solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
  theta <- solvebtb%*% ttt
  cat("Diffusion tensors generated \n")

  res <- ttt - t(from@btb) %*% theta
  mres2 <- res[1,]^2
  for(i in 2:ngrad) mres2 <- mres2 + res[i,]^2
  sigma2 <- array(mres2/(ngrad-6),ddim)
  dim(theta) <- c(6,ddim)
  dim(res) <- c(ngrad,ddim)
  cat("Variance estimates generated \n")

  rm(mres2)
  gc()

  dim(s0) <- ddim
  mask <- s0>from@level
  lags <- c(5,5,3)
  scorr <- .Fortran("mcorr",as.double(aperm(res,c(2:4,1))),
                   as.logical(mask),
                   as.integer(ddim[1]),
                   as.integer(ddim[2]),
                   as.integer(ddim[3]),
                   as.integer(ngrad),
                   scorr = double(prod(lags)),
                   as.integer(lags[1]),
                   as.integer(lags[2]),
                   as.integer(lags[3]),
                   PACKAGE="dti",DUP=FALSE)$scorr
  dim(scorr) <- lags
  scorr[is.na(scorr)] <- 0
  cat("first order  correlation in x-direction",signif(scorr[2,1,1],3),"\n")
  cat("first order  correlation in y-direction",signif(scorr[1,2,1],3),"\n")
  cat("first order  correlation in z-direction",signif(scorr[1,1,2],3),"\n")
  bw <- optim(c(2,2,2),corrrisk,method="L-BFGS-B",lower=c(.25,.25,.25),lag=lags,data=scorr)$par
  bw[bw <= .25] <- 0

  invisible(new(to,
                list(theta = theta, sigma = sigma2, scorr = scorr, bw = bw),
                btb   = from@btb,
                ngrad = from@ngrad, # = dim(btb)[2]
                ddim  = from@ddim,
                ddim0 = from@ddim0,
                xind  = from@xind,
                yind  = from@yind,
                zind  = from@zind,
                source= from@source)
            )
})

dtiTensor <- function(object, ...) cat("No DTI tensor calculation defined for this class:",class(object),"\n")

setGeneric("dtiTensor", function(object, ...) standardGeneric("dtiTensor"))

setMethod("dtiTensor","dtiData",
function(object) {
  ngrad <- object@ngrad
  ddim <- object@ddim

  s0 <- object$s0
  si <- object$si
  dim(s0) <- dim(si) <- NULL
  ttt <- -log(si/s0)
  ttt[is.na(ttt)] <- 0
  ttt[(ttt == Inf)] <- 0
  ttt[(ttt == -Inf)] <- 0
  dim(ttt) <- c(prod(ddim),ngrad)
  ttt <- t(ttt)
  cat("Data transformation completed \n")

  btbsvd <- svd(object@btb)
  solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
  theta <- solvebtb%*% ttt
  cat("Diffusion tensors generated \n")

  res <- ttt - t(object@btb) %*% theta
  mres2 <- res[1,]^2
  for(i in 2:ngrad) mres2 <- mres2 + res[i,]^2
  sigma2 <- array(mres2/(ngrad-6),ddim)
  dim(theta) <- c(6,ddim)
  dim(res) <- c(ngrad,ddim)
  cat("Variance estimates generated \n")

  rm(mres2)
  gc()

  dim(s0) <- ddim
  mask <- s0 > object@level
  lags <- c(5,5,3)
  scorr <- .Fortran("mcorr",as.double(aperm(res,c(2:4,1))),
                   as.logical(mask),
                   as.integer(ddim[1]),
                   as.integer(ddim[2]),
                   as.integer(ddim[3]),
                   as.integer(ngrad),
                   scorr = double(prod(lags)),
                   as.integer(lags[1]),
                   as.integer(lags[2]),
                   as.integer(lags[3]),
                   PACKAGE="dti",DUP=FALSE)$scorr
  dim(scorr) <- lags
  scorr[is.na(scorr)] <- 0
  cat("first order  correlation in x-direction",signif(scorr[2,1,1],3),"\n")
  cat("first order  correlation in y-direction",signif(scorr[1,2,1],3),"\n")
  cat("first order  correlation in z-direction",signif(scorr[1,1,2],3),"\n")

  scorr[is.na(scorr)] <- 0
  bw <- optim(c(2,2,2),corrrisk,method="L-BFGS-B",lower=c(.25,.25,.25),lag=lags,data=scorr)$par
  bw[bw <= .25] <- 0

  invisible(new("dtiTensor",
                list(theta = theta, sigma = sigma2, scorr = scorr, bw = bw),
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                ddim  = object@ddim,
                ddim0 = object@ddim0,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                source= object@source)
            )
})

create.designmatrix.dti <- function(gradient, bvalue=1) {
  dgrad <- dim(gradient)
  if (dgrad[2]==3) gradient <- t(gradient)
  if (dgrad[1]!=3) stop("Not a valid gradient matrix")

  btb <- matrix(0,6,dgrad[2])
  btb[1,] <- gradient[1,]*gradient[1,]
  btb[4,] <- gradient[2,]*gradient[2,]
  btb[6,] <- gradient[3,]*gradient[3,]
  btb[2,] <- 2*gradient[1,]*gradient[2,]
  btb[3,] <- 2*gradient[1,]*gradient[3,]
  btb[5,] <- 2*gradient[2,]*gradient[3,]

  btb * bvalue
}

# really setAs() or setMethod?
setAs("dtiTensor","dtiIndices",function(from,to) {
  ddim <- from@ddim

  ll <- array(0,c(ddim,3))
  th <- array(0,c(ddim,9))
  ierr <- array(0,ddim)

  for (i in 1:ddim[1]) {
    cat(".")
    for (j in 1:ddim[2]) {
      for (k in 1:ddim[3]) {
        z <- .Fortran("eigen3",
                      as.double(from$theta[,i,j,k]),
                      lambda = double(3),
                      theta = double(3*3),
                      ierr = integer(1),
                      PACKAGE="dti")[c("lambda","theta","ierr")]
        ll[i,j,k,] <- z$lambda
        th[i,j,k,] <- z$theta
        ierr[i,j,k] <- z$ierr
      }
    }
  }
  cat("\ncalculated eigenvalues and -vectors\n")

  dim(th) <- c(ddim,3,3)
  dim(ll) <- c(prod(ddim),3)

  trc <- as.vector(ll %*% c(1,1,1))/3
  fa <- sqrt(1.5*((sweep(ll,1,trc)^2)%*% c(1,1,1))/((ll^2)%*% c(1,1,1)))
  ra <- sqrt(((sweep(ll,1,trc)^2)%*% c(1,1,1))/(3*trc))

  cat("calculated anisotropy indices\n")

  bary <- c((ll[,3] - ll[,2]) / (3*trc) , 2*(ll[,2] - ll[,1]) / (3*trc) , ll[,1] / trc)

  dim(ll) <- c(ddim,3)
  dim(trc) <- dim(fa) <- dim(ra) <- ddim
  dim(bary) <- c(ddim,3)

  cat("calculated barycentric coordinates\n")

  invisible(new(to,
                list(fa = fa, ra = ra, trc = trc, bary = bary, lambda = ll, eigenv = th),
                btb   = from@btb,
                ngrad = from@ngrad, # = dim(btb)[2]
                ddim  = from@ddim,
                ddim0 = from@ddim0,
                xind  = from@xind,
                yind  = from@yind,
                zind  = from@zind,
                source= from@source)
            )


})

dtiIndices <- function(object, ...) cat("No DTI indices calculation defined for this class:",class(object),"\n")

setGeneric("dtiIndices", function(object, ...) standardGeneric("dtiIndices"))

setMethod("dtiIndices","dtiTensor",
function(object, which) {
  ddim <- object@ddim

  ll <- array(0,c(ddim,3))
  th <- array(0,c(ddim,9))
  ierr <- array(0,ddim)

  for (i in 1:ddim[1]) {
    cat(".")
    for (j in 1:ddim[2]) {
      for (k in 1:ddim[3]) {
        z <- .Fortran("eigen3",
                      as.double(object$theta[,i,j,k]),
                      lambda = double(3),
                      theta = double(3*3),
                      ierr = integer(1),
                      PACKAGE="dti")[c("lambda","theta","ierr")]
        ll[i,j,k,] <- z$lambda
        th[i,j,k,] <- z$theta
        ierr[i,j,k] <- z$ierr
      }
    }
  }
  cat("\ncalculated eigenvalues and -vectors\n")

  dim(th) <- c(ddim,3,3)
  dim(ll) <- c(prod(ddim),3)

  trc <- as.vector(ll %*% c(1,1,1))/3
  cat("voxel with negative trace",sum(trc<=0),"\n")
  ind <- trc > 0
  fa <- ra <- numeric(prod(ddim))
  fa[ind] <- sqrt((1.5*((sweep(ll,1,trc)^2)%*% c(1,1,1))/((ll^2)%*% c(1,1,1)))[ind])
  ra[ind] <- sqrt(((sweep(ll,1,trc)^2)%*% c(1,1,1))[ind]/(3*trc[ind]))
# set fa, ra to zero if trc <= 0
  cat("calculated anisotropy indices\n")

  bary <- c((ll[,3] - ll[,2]) / (3*trc) , 2*(ll[,2] - ll[,1]) / (3*trc) , ll[,1] / trc)

  dim(ll) <- c(ddim,3)
  dim(trc) <- dim(fa) <- dim(ra) <- ddim
  dim(bary) <- c(ddim,3)

  cat("calculated barycentric coordinates\n")

  invisible(new("dtiIndices",
                list(fa = fa, ra = ra, trc = trc, bary = bary, lambda = ll, eigenv = th),
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                ddim  = object@ddim,
                ddim0 = object@ddim0,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                source= object@source)
            )


})

