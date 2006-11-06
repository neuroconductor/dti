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
     title(paste("Anisotropy index (h=",signif(hakt,3),"), slice",slice,"range:",signif(range(z$anindex[z$mask]))))
     image(z$bi[,,slice],col=grey((0:255)/255))
     title(paste("sum of weights  max=",signif(max(z$bi),3),"mean=",signif(mean(z$bi[z$mask]),3)))
     }
     hakt <- hakt*hincr
  }
invisible(list(theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,mask=z$mask,call=args))
}
dtianiso2<-function(y,hmax,lambda,rho,graph=FALSE,slice=NULL,bvec=NULL,sigma2=NULL,scorr=c(.5,.5),mask=NULL,quant=.8,zext=1){
  args <- match.call()
  btb<-matrix(0,6,dim(bvec)[2])
  btb[1,]<-bvec[1,]^2
  btb[4,]<-bvec[2,]^2
  btb[6,]<-bvec[3,]^2
  btb[2,]<-2*bvec[1,]*bvec[2,]
  btb[3,]<-2*bvec[1,]*bvec[3,]
  btb[5,]<-2*bvec[3,]*bvec[2,]
  Bcov <- btb%*%t(btb)
  dimy <- dim(y)
  if(length(dimy)!=4||dimy[1]!=6) stop("y does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  if(is.null(mask)) mask <- array(logical(n),dimy[-1])
  if(is.null(dim(sigma2))) {
    sigma2 <- rep(sigma2,n)
    dim(sigma2) <- dimy[-1]
  }
  sigma2[sigma2<=mean(sigma2)*1e-5]<- mean(sigma2)*1e-5
  z <- .Fortran("projdt",
                as.double(y),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                ynew=double(6*n),
                mask=logical(n),
                DUP=FALSE,
                PACKAGE="dti")[c("ynew","mask")]
  y <- array(z$ynew,dimy)
  mask <- array(z$mask,dimy[-1])&mask
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
  z$bi <- 1/sigma2
  dim(z$theta) <- dimy
  dim(z$anindex) <-dim(z$det) <-dim(z$mask) <- dimy[-1]
  z$mask <- array(z$mask,dimy[-1])&mask
  dim(z$andirection) <- c(3,dimy[-1]) 
#  initial state for h=1
#
  if(graph){
     require(adimpro)
     oldpar <- par(mfrow=c(1,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
     on.exit(par(oldpar))
     if(is.null(slice)) slice<-n3%/%2
     show.image(make.image(andir.image(z,quant=quant)[,,slice,]))
     title(paste("Anisotropy index (h=1), slice",slice))
     ni<-z$bi[,,slice]*sigma2[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights, slice",slice))
  }
  if (max(scorr)>0) {
    h0 <- numeric(length(scorr))
    for (i in 1:length(h0)) h0[i] <- get.bw.gauss(scorr[i],interv=2)
    if (length(h0)<2) h0 <- rep(h0[1],2)
# no spatial correlation iz z-direction
    cat("Corresponding bandwiths for specified correlation:",h0,"\n")
  }
  hincr <- 1.25^(1/3)
  hakt0 <- 1
  hakt <- hincr
  lambda0 <- lambda
  while( hakt <= hmax) {
    if (scorr[1]>=0.1) lambda0 <- lambda * Spatialvar.gauss(hakt0/0.42445/4,h0,2) /
      Spatialvar.gauss(h0,1e-5,2) /
        Spatialvar.gauss(hakt0/0.42445/4,1e-5,2)
     z <- .Fortran("awsdti2",
                as.double(y),
                as.double(z$theta),
                bi=as.double(z$bi),
                anindex=as.double(z$anindex),
                andirection=as.double(z$andirection),
                det=as.double(z$det),
                as.double(Bcov),
                as.double(sigma2),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.double(hakt),
                as.double(zext),
                as.double(rho),
                as.double(lambda0),
                theta=double(6*n),
                mask=as.logical(z$mask),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","bi","anindex","andirection","det","mask")]
     dim(z$bi) <- dim(z$anindex) <-dim(z$det) <- dimy[-1]
     dim(z$theta) <- dimy
     dim(z$andirection) <- c(3,dimy[-1]) 
     if(graph){
     show.image(make.image(andir.image(z,quant=quant)[,,slice,]))
     title(paste("Anisotropy index (h=",signif(hakt,3),"), slice",slice,"range:",signif(min(z$anindex[z$mask]),3),"-",
                                                                                 signif(max(z$anindex[z$mask]),3)))
     ni<-z$bi[,,slice]*sigma2[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights  mean=",signif(mean(z$bi[z$mask]*sigma2[z$mask]),3)))
     }
     cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex[z$mask],c(.5, .75, .9, .95, 1)),3),"\n")
     hakt0<-hakt
     hakt <- hakt*hincr
  }
invisible(list(theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,mask=z$mask,call=args))
}

andir.image <- function(dtobject,method=1,quant=0){
anindex <- dtobject$anindex
andirection <- dtobject$andirection
mask <- dtobject$mask
anindex[anindex>1]<-0
anindex[anindex<0]<-0
dimg <- dim(anindex)
dim(andirection)<-c(3,prod(dimg))
minanindex <- quantile(anindex[mask],quant)
if(method==1) {
andirection[1,] <- abs(andirection[1,])
andirection[2,] <- abs(andirection[2,])
andirection[3,] <- abs(andirection[3,])
} else {
ind<-andirection[1,]<0
andirection[,ind] <- - andirection[,ind]
andirection[2,] <- (1+andirection[2,])/2
andirection[3,] <- (1+andirection[3,])/2
}
andirection <- t(andirection)
andirection <- andirection*as.vector(anindex)*as.numeric(mask)*as.numeric(anindex>minanindex)
dim(andirection)<-c(dimg,3)
invisible(andirection)
} 

Spatialvar.gauss<-function(h,h0,d,interv=1){
#
#   Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h in 
#
#   case of colored noise that was produced by smoothing with Gaussian kernel and bandwidth h0
#
#   Spatialvar.gauss(lkern,h,h0,d)/Spatialvar.gauss(lkern,h,1e-5,d) gives the 
#   a factor for lambda to be used with bandwidth h 
#
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing, 
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding 
#  discretisation into voxel) 
#
  h0 <- pmax(h0,1e-5)
  h <- pmax(h,1e-5)
  h<-h/2.3548*interv
  if(length(h)==1) h<-rep(h,d)
  ih<-trunc(4*h)
  ih<-pmax(1,ih)
  dx<-2*ih+1
  penl<-dnorm(((-ih[1]):ih[1])/h[1])
  if(d==2) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),dnorm(((-ih[2]):ih[2])/h[2]),"*")
  if(d==3) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
  dim(penl)<-dx
  h0<-h0/2.3548*interv
  if(length(h0)==1) h0<-rep(h0,d)
  ih<-trunc(4*h0)
  ih<-pmax(1,ih)
  dx0<-2*ih+1
  x<- ((-ih[1]):ih[1])/h0[1]
  penl0<-dnorm(((-ih[1]):ih[1])/h0[1])
  if(d==2) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),dnorm(((-ih[2]):ih[2])/h0[2]),"*")
  if(d==3) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),outer(dnorm(((-ih[2]):ih[2])/h0[2]),dnorm(((-ih[3]):ih[3])/h0[3]),"*"),"*")
  dim(penl0)<-dx0
  penl0<-penl0/sum(penl0)
  dz<-dx+dx0-1
  z<-array(0,dz)
  if(d==1){
    for(i1 in 1:dx0) {
      ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
      ind1<-ind1[ind1<=dz][-1]
      z[-ind1]<-z[-ind1]+penl*penl0[i1]
    }
  } else if(d==2){
    for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]){
      ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
      ind1<-ind1[ind1<=dz[1]][-1]
      ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
      ind2<-ind2[ind2<=dz[2]][-1]
      z[-ind1,-ind2]<-z[-ind1,-ind2]+penl*penl0[i1,i2]
    }
  } else if(d==3){
    for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]) for(i3 in 1:dx0[3]){
      ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
      ind1<-ind1[ind1<=dz[1]][-1]
      ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
      ind2<-ind2[ind2<=dz[2]][-1]
      ind3<-c(0:(i3-1),(dz[3]-dx0[3]+i3):dz[3]+1)
      ind3<-ind3[ind3<=dz[3]][-1]
      z[-ind1,-ind2,-ind3]<-z[-ind1,-ind2,-ind3]+penl*penl0[i1,i2,i3]
    }
  }
  sum(z^2)/sum(z)^2*interv^d
}
get.bw.gauss <- function(corr, step = 1.001,interv=2) {
  
  # get the   bandwidth for lkern corresponding to a given correlation
  #  keep it simple result does not depend on d

  #  interv allows for further discretization of the Gaussian Kernel, result depends on
  #  interv for small bandwidths. interv=1  is correct for kernel smoothing, 
  #  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding 
  #  discretisation into voxel)   
  if (corr < 0.1) {
    h <- 0
  } else { 
    h <- .5
    z <- 0
    while (z<corr) {
      h <- h*step
      z <- get.corr.gauss(h,interv)
    }
    h <- h/step
  }
  h
}

get.corr.gauss <- function(h,interv=1) {
    #
    #   Calculates the correlation of 
    #   colored noise that was produced by smoothing with "gaussian" kernel and bandwidth h
    #   Result does not depend on d for "Gaussian" kernel !!
    h <- h/2.3548*interv
    ih <- trunc(4*h+ 2*interv-1)
    dx <- 2*ih+1
    penl <- dnorm(((-ih):ih)/h)
    sum(penl[-(1:interv)]*penl[-((dx-interv+1):dx)])/sum(penl^2)
}
