### CODE IN THIS FILE IS STILL EXPERIMENTAL ###

create.designmatrix.dti <- function(bvec, bvalue=1) {
  cat("\nNOTE: This code is still experimental!\n") 
  dimension <- dim(bvec)[2] # should be 3
  if (dimension != 3) {
    warning("Error: gradient vectors do not have length 3")
    return(invisible(NULL))
  }
  directions <- dim(bvec)[1] # number of measured directions

  z <- matrix(0, directions, 6)

  for (d in 1:directions) {
    z[d,1] <- bvec[d,1]^2
    z[d,4] <- bvec[d,2]^2
    z[d,6] <- bvec[d,3]^2
    z[d,2] <- 2*bvec[d,1]*bvec[d,2]
    z[d,3] <- 2*bvec[d,1]*bvec[d,3]
    z[d,5] <- 2*bvec[d,2]*bvec[d,3]
    z[d,] <- bvalue*z[d,]
  }
  
  z
}

calculate.lm.dti <- function(ttt,z,res=FALSE) {
  cat("\nNOTE: This code is still experimental!\n") 
  svdresult <- svd(z)
  u <- svdresult$u
  v <- svdresult$v
  vt <- t(v)
  lambda1 <- diag(1/svdresult$d)
  lambda2 <- diag(1/svdresult$d^2)
  xtx <- v %*% lambda2 %*% vt
  # now we have z = u lambda1^(-1) vt
  
  # define some variables and make ttt a matrix
  dy <- dim(ttt)
  voxelcount <- prod(dy[1:3])
  dim(ttt) <- c(prod(dy[1:3]),dy[4])
  
  # calculate the paramters and residuals for all voxels
  beta <- ttt %*% u %*% lambda1 %*% vt
  residuals <- ttt - beta %*% t(z)
  b <- rep(1/dy[4],length=dy[4])
  variance <- ((residuals^2 %*% b) * dim(z)[1] / (dim(z)[1]-dim(z)[2]))
  dim(variance) <- c(dy[1:3])
  dim(beta) <- c(dy[1:3],dim(z)[2])
  dim(residuals) <- c(dy[1:3],dim(z)[1])

  if (res) {
    list(dt=beta,xtx=xtx,variance=variance,residuals=residuals)
  } else {
    list(dt=beta,xtx=xtx,variance=variance)    
  }
}
  
determine.eigenvalue <- function(y, reduced=FALSE) {
  cat("\nNOTE: This code is still experimental!\n") 

  dy <- dim(y)
  n <- prod(dy[-4]) # number of voxel

  if (reduced) {
    ll <- array(0,dy[1:3])
    th <- array(0,c(dy[1:3],3))
  } else {
    ll <- array(0,c(dy[1:3],3))
    th <- array(0,c(dy[1:3],9))
  }
  ierr <- array(0,dy[1:3])

  for (i in 1:dy[1]) {
    cat(".")
    for (j in 1:dy[2]) {
      for (k in 1:dy[3]) {
        if (reduced) {
          z <- .Fortran("eigen3r",
                        as.double(y[i,j,k,]),
                        lambda = double(1),
                        theta = double(3),
                        ierr = integer(1),
                        PACKAGE="dti")[c("lambda","theta","ierr")]
          ll[i,j,k] <- z$lambda
        } else {
          z <- .Fortran("eigen3",
                        as.double(y[i,j,k,]),
                        lambda = double(3),
                        theta = double(3*3),
                        ierr = integer(1),
                        PACKAGE="dti")[c("lambda","theta","ierr")]
          ll[i,j,k,] <- z$lambda
        }
        th[i,j,k,] <- z$theta
        ierr[i,j,k] <- z$ierr
      }
    }
  }

  if (!reduced) dim(th) <- c(dy[1:3],3,3)
  
  list(lambda = ll, theta = th, ierr = ierr)
}

anisotropy <- function(eigen) {
  cat("\nNOTE: This code is still experimental!\n") 
  dimdt <- dim(eigen)
  dim(eigen) <- c(prod(dimdt[1:3]),3)

  trc <- as.vector(eigen %*% c(1,1,1))/3
  fa <- sqrt(1.5*((sweep(eigen,1,trc)^2)%*% c(1,1,1))/((eigen^2)%*% c(1,1,1)))
  ra <- sqrt(((sweep(eigen,1,trc)^2)%*% c(1,1,1))/(3*trc))

  dim(trc) <- dim(fa) <- dim(ra) <- dimdt[1:3]
  
  list(fa=fa, ra=ra, trace=trc)
}

theta.estimate <- function(y,dt=NULL,h) {
  cat("\nNOTE: This code is still experimental!\n")

  n1 <- dim(y)[1]
  n2 <- dim(y)[2]
  n3 <- dim(y)[3]

  if (is.null(dt)) dt <- array(1,c(n1,n2,n3))
  
  z <- .Fortran("esttheta",
                as.double(aperm(y,c(4,1,2,3))),
                as.double(dt),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.double(h),
                theta = double(3*n1*n2*n3),
                PACKAGE="dti")[c("theta")]

  dim(z$theta) <- c(3,n1,n2,n3)
  z$theta
}

dt.estimate <- function(theta,y,hd,ht) {
  cat("\nNOTE: This code is still experimental!\n") 

  try <- dftr(y)
  cat("Trace determined\n")
  
  n1 <- dim(try)[1]
  n2 <- dim(try)[2]
  n3 <- dim(try)[3]
  
  if (length(theta) == 3) theta <- aperm(array(theta,dim=c(3,n1,n2,n3)),c(2,3,4,1))
  cat("Extended theta\n")
  
  z <- .Fortran("estimdt",
                as.double(theta),
                as.double(try),
                as.double(hd),
                as.double(ht),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                dt = double(n1*n2*n3),
                PACKAGE="dti")[c("dt")]

  dim(z$dt) <- c(n1,n2,n3)
  z$dt
}

dftr <- function(y) {
  dy <- dim(y)
  dim(y) <- c(prod(dy[1:3]),6)
  try <- y[,c(1,4,6)]%*%c(1,1,1)
  dim(try) <- dy[1:3]
  try
}

dt.estimate2 <- function(theta,y,ai,lmax,hd,ht) {
  cat("\nNOTE: This code is still experimental!\n") 

  try <- dftr(y)
  cat("Trace determined\n")
  
  n1 <- dim(try)[1]
  n2 <- dim(try)[2]
  n3 <- dim(try)[3]
  
  if (length(theta) == 3) theta <- aperm(array(theta,dim=c(3,n1,n2,n3)),c(2,3,4,1))
  cat("Extended theta\n")
  
  z <- .Fortran("estimdt2",
                as.double(theta),
                as.double(try),
		as.double(ai),
		as.double(lmax),
                as.double(hd),
                as.double(ht),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                dt = double(n1*n2*n3),
                PACKAGE="dti")[c("dt")]

  dim(z$dt) <- c(n1,n2,n3)
  z$dt
}

tensor.estimate <- function(y,dt=NULL,h) {
  cat("\nNOTE: This code is still experimental!\n")

  n1 <- dim(y)[1]
  n2 <- dim(y)[2]
  n3 <- dim(y)[3]
  n <- n1*n2*n3
  if (is.null(dt)) dt <- array(1,c(n1,n2,n3))
  
  z <- .Fortran("esttens",
                as.double(aperm(y,c(4,1,2,3))),
                as.double(dt),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.double(h),
                theta = double(3*n),
		yhat = double(6*n),
		aihat = double(n),
		lmhat = double(n),
                PACKAGE="dti")[c("theta","yhat","aihat","lmhat")]
              
  dim(z$theta) <- c(3,n1,n2,n3)
  dim(z$yhat) <- c(6,n1,n2,n3)
  z$yhat <- aperm(z$yhat,c(2:4,1))
  dim(z$aihat) <- c(n1,n2,n3)
  dim(z$lmhat) <- c(n1,n2,n3)
  z
}

dtianiso<-function(y,hmax,lambda,rho,graph=FALSE,slice=NULL){
  args <- match.call()
  dimy <- dim(y)
  if(length(dimy)!=4||dimy[1]!=6) stop("y does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  theta <- y
  z <- .Fortran("initdti",
                theta=as.double(theta),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                mask=logical(n),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","anindex","andirection","det","mask")]
  z$bi <- array(1,dimy[-1])
  dim(z$theta) <- dimy
  dim(z$anindex) <-dim(z$det) <-dim(z$mask) <- dimy[-1]
  dim(z$andirection) <- c(3,dimy[-1])
#
#  initial state for h=1
#
  if(graph){
     oldpar <- par(mfrow=c(1,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
     on.exit(par(oldpar))
     if(is.null(slice)) slice<-n3%/%2
     image(z$anindex[,,slice],col=grey((0:255)/255))
     title(paste("Anisotropy index (h=1), slice",slice))
     image(z$bi[,,slice],col=grey((0:255)/255))
     title(paste("sum of weights, slice",slice))
  }
  hincr <- 1.25^(1/3)
  hakt <- hincr
  while( hakt <= hmax) {
     z <- .Fortran("awsdti",
                as.double(y),
                as.double(z$theta),
                bi=as.double(z$bi),
                anindex=as.double(z$anindex),
                andirection=as.double(z$andirection),
                det=as.double(z$det),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.double(hakt),
                as.double(rho),
                as.double(lambda),
                theta=double(6*n),
                mask=as.logical(z$mask),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","bi","anindex","andirection","det","mask")]
     dim(z$bi) <- dim(z$anindex) <-dim(z$det) <- dimy[-1]
     dim(z$theta) <- dimy
     dim(z$andirection) <- c(3,dimy[-1])
     if(graph){
     image(z$anindex[,,slice],col=grey((0:255)/255))
     title(paste("Anisotropy index (h=",signif(hakt,3),"), slice",slice))
     image(z$bi[,,slice],col=grey((0:255)/255))
     title(paste("sum of weights  max=",signif(max(z$bi),3),"mean=",signif(mean(z$bi[z$mask]),3)))
     }
     hakt <- hakt*hincr
  }
invisible(list(theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,call=args))
}