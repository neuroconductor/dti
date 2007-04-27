  # solvebtb mit einfuegen??

setClass("dti",
         representation(btb    = "matrix",
                        ngrad  = "integer", # = dim(btb)[2]
                        ddim   = "integer",
                        ddim0  = "integer",
                        xind   = "integer",
                        yind   = "integer",
                        zind   = "integer",
                        source = "character")
         )

setClass("dtiData",
         representation(level = "numeric"),
         contains=c("list","dti"),
         validity=function(object){
          if (sum(c("s0","si") %in% names(object)) != 2) {
            cat("s0 and/or si not in data\n")
            return(invisible(FALSE))
          }
          if (!("array" %in% class(object$s0))) {
            cat("s0 not an array\n")
            return(invisible(FALSE))
          }
          if (!("array" %in% class(object$si))) {
            cat("s0 not an array\n")
            return(invisible(FALSE))
          }
          if (length(dim(object$s0)) != 3) {
            cat("dimension of s0 is",dim(object$s0),", but we want 3 dimensions\n")
            return(invisible(FALSE))
          }
          if (length(dim(object$si)) != 4) {
            cat("dimension of si is",dim(object$si),", but we want 4 dimensions\n")
            return(invisible(FALSE))
          }
          if (dim(object$s0) != object@ddim) {
            cat("dimension of s0 is",dim(object$s0),", but we want",object@ddim,"\n")
            return(invisible(FALSE))
          }
          if (!all(dim(object$si) != c(object@ddim,object@ngrad))) {
            cat("dimension of si is",dim(object$si),", but we want",c(object@ddim,object@ngrad),"\n")
            return(invisible(FALSE))
          }
         }
         )
setClass("dtiTensor",
         contains=c("list","dti"),
         validity=function(object){
          if (sum(c("theta","sigma","scorr") %in% names(object)) != 3) {
            cat("s0 and/or si not in data\n")
            return(invisible(FALSE))
          }
          if (!("array" %in% class(object$theta))) {
            cat("s0 not an array\n")
            return(invisible(FALSE))
          }
          if (!("array" %in% class(object$sigma))) {
            cat("s0 not an array\n")
            return(invisible(FALSE))
          }
         }
         )

setClass("dtiIndices",
         representation(fa     = "array",
                        ra     = "array",
                        trc    = "array",
                        bary   = "array",
                        lambda = "array",
                        eigenv = "array"),
         contains=c("list","dti"),
          validity=function(object){
          if (sum(c("fa","ra","trc","lambda","eigenv") %in% names(object)) != 5) {
            cat("fa,ra, trc, lambda, or eigenv not in data\n")
            return(invisible(FALSE))
          }
         }
        )




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
function(x, y, slice=1, method=1, quant=0, minanindex=NULL, show=TRUE, ...) {
  cat("Plot called with class",class(x),"\n")
  if(is.null(x$fa)) cat("No anisotropy index yet")
  adimpro <- require(adimpro)
  anindex <- x$fa[,,slice]
  dimg <- x@ddim[1:2]
cat("A\n")
  andirection <- x$eigenv[,,slice,,]
  anindex[anindex>1]<-0
  anindex[anindex<0]<-0
  dim(andirection)<-c(prod(dimg),3,3)
  if(is.null(minanindex)) minanindex <- quantile(anindex,quant,na.rm=TRUE)
  if(method==1) {
    andirection[,1,3] <- abs(andirection[,1,3])
    andirection[,2,3] <- abs(andirection[,2,3])
    andirection[,3,3] <- abs(andirection[,3,3])
  } else {
    ind<-andirection[,1,3]<0
    andirection[ind,,] <- - andirection[ind,,]
    andirection[,2,3] <- (1+andirection[,2,3])/2
    andirection[,3,3] <- (1+andirection[,3,3])/2
  }
cat("B",dim(andirection),"\n")
  andirection <- andirection[,,1]
cat("C",dim(andirection),"\n")
cat("D",length(anindex),"\n")
  andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
cat("E",dim(andirection),"\n")
cat("F",dimg,"\n")
  dim(andirection)<-c(dimg,3)
  if(adimpro) {
cat("G",range(andirection),"\n")
cat("H",range(andirection,na.rm=TRUE),"\n")
andirection[is.na(andirection)] <- 0
    andirection <- make.image(andirection)
    if(show) show.image(andirection,...)
  } else if(show) {
    dim(anindex) <- dimg
    image(anindex,...)
  }
  invisible(andirection)
})


dtiData <- function(gradient,imagefile,ddim,xind=NULL,yind=NULL,zind=NULL,level=0) {
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
  scorr <- c(0,0)
  res <- aperm(res,c(2:4,1))
  dim(res) <- c(prod(ddim),ngrad)
  res1 <- as.vector(res[as.vector(mask),])
  scorr[1] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
  cat("correlation in x-direction",signif(scorr[1],3),"\n")
  dim(res) <- c(ddim,ngrad)
  res <- aperm(res,c(2,1,3,4))
  dim(res) <- c(prod(ddim),ngrad)
  res1 <- as.vector(res[as.vector(aperm(mask,c(2,1,3))),])
  scorr[2] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
  cat("correlation in y-direction",signif(scorr[2],3),"\n")


  invisible(new(to,
                list(theta = theta, sigma = sigma2, scorr = scorr),
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
  mask <- s0>object@level
  scorr <- c(0,0)
  res <- aperm(res,c(2:4,1))
  dim(res) <- c(prod(ddim),ngrad)
  res1 <- as.vector(res[as.vector(mask),])
  scorr[1] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
  cat("correlation in x-direction",signif(scorr[1],3),"\n")
  dim(res) <- c(ddim,ngrad)
  res <- aperm(res,c(2,1,3,4))
  dim(res) <- c(prod(ddim),ngrad)
  res1 <- as.vector(res[as.vector(aperm(mask,c(2,1,3))),])
  scorr[2] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
  cat("correlation in y-direction",signif(scorr[2],3),"\n")


  invisible(new("dtiTensor",
                list(theta = theta, sigma = sigma2, scorr = scorr),
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

  bary <- c((ll[,1] - ll[,2]) / trc , 2*(ll[,2] - ll[,3]) / trc , 3*ll[,3] / trc)

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
  fa <- sqrt(1.5*((sweep(ll,1,trc)^2)%*% c(1,1,1))/((ll^2)%*% c(1,1,1)))
  ra <- sqrt(((sweep(ll,1,trc)^2)%*% c(1,1,1))/(3*trc))

  cat("calculated anisotropy indices\n")

  bary <- c((ll[,1] - ll[,2]) / trc , 2*(ll[,2] - ll[,3]) / trc , 3*ll[,3] / trc)

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

