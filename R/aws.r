# This file contains the implementation of dti.smooth() for 
# "dtiData" lines 12--
# and 
# "dtiTensor" lines 225--
# objects. The latter can be done with 
# "Euclidian Metric" lines 237--
# "Riemann Metric" lines 414--

dti.smooth <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dti.smooth", function(object, ...) standardGeneric("dti.smooth"))

setMethod("dti.smooth", "dtiData", function(object,hmax=5,hinit=NULL,lambda=25,
                                            rho=1,graph=FALSE,slice=NULL,quant=.8,
                                            minanindex=NULL,zext=1,eps=1e-6,hsig=2.5) {
  if (graph) {
    adimpro <- require(adimpro)
    if (!adimpro) cat("No graphical output! Install package adimpro from CRAN!\n")
    graph <- graph & adimpro
  }
  args <- match.call()
  si <- object$si
  s0 <- object$s0
  ngrad <- object@ngrad
  ddim0 <- object@ddim0
  ddim <- object@ddim
  xind <- object@xind
  yind <- object@yind
  zind <- object@zind
  source <- object@source
  btb <- object@btb
  Bcov <- btb%*%t(btb)
  btbsvd <- svd(btb)
  solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
  projectmat <- diag(ngrad) - btbsvd$v %*% t(btbsvd$v)

  dtobject <- as(object,"dtiTensor")
  y <- dtobject$theta
  sigma2 <- dtobject$sigma
  scorr <- dtobject$scorr
  h0 <- dtobject$bw
  cat("Corresponding bandwiths for specified correlation:",h0,"\n")

  rm(object,dtobject)
  gc()
  dimy <- dim(y)
  if(length(dimy)!=4||dimy[1]!=6) stop("y does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  sigma2[sigma2<=mean(sigma2)*1e-5]<- mean(sigma2)*1e-5
  z <- .Fortran("projdt2",
                as.double(y),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                theta=double(6*n),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                as.double(eps),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","anindex","andirection","det")]
  y <- array(z$theta,dimy)
  z$bi <- 1/sigma2
  dim(z$theta) <- dimy
  dim(z$anindex) <-dim(z$det) <- dimy[-1]
  dim(z$andirection) <- c(3,dimy[-1]) 
  z$s0hat <- s0
  z <- .Fortran("smsigma",
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.double(hsig),
                       as.double(zext),
                       sigma2hat=double(n1*n2*n3),
                       DUP=FALSE,
                       PACKAGE="dti")["sigma2hat"]
#
#  initial state for h=1
#
  if(graph){
     oldpar <- par(mfrow=c(3,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
     on.exit(par(oldpar))
     if(is.null(slice)) slice<-n3%/%2
     class(z) <- "dti"
     img<-z$theta[1,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dxx: mean",signif(mean(img),3),"max",signif(max(img),3)))
     img<-z$theta[2,,,slice]
     show.image(make.image(img))
     title(paste("Dxy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[3,,,slice]
     show.image(make.image(img))
     title(paste("Dxz: min",signif(min(img),3),"max",signif(max(img),3)))
     show.image(make.image(z$anindex[,,slice]))
     title(paste("Anisotropy index  range:",signif(min(z$anindex),3),"-",
                  signif(max(z$anindex),3)))
     img<-z$theta[4,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: mean",signif(mean(img),3),"max",signif(max(img),3)))
     img<-z$theta[5,,,slice]
     show.image(make.image(img))
     title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
     plot(c(0,1),c(0,1))
     title(paste("Directions (h=1), slice",slice))
     ni<-z$bi[,,slice]*sigma2[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights  mean=",signif(mean(z$bi*sigma2),3)))
     img<-z$theta[6,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dzz: mean",signif(mean(img),3),"max",signif(max(img),3)))
  }
  hincr <- 1.25^(1/3)
  if(is.null(hinit)){
  hakt0 <- 1
  hakt <- hincr
  } else {
  hakt0 <- max(1,hinit/hincr)
  hakt <- hinit
  }
  lambda0 <- lambda
  while(hakt <= hmax) {
    if (any(h0 >= 0.25)) {
       corrfactor <- Spatialvar.gauss(hakt0/0.42445/4,h0,3) /
       Spatialvar.gauss(h0,1e-5,3) /
       Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
       lambda0 <- lambda * corrfactor
       cat("Correction factor for spatial correlation",signif(corrfactor,3),"\n")
    }
     z <- .Fortran("awssidti",
                as.double(s0),
                as.double(si),
                as.double(z$theta),
                bi=as.double(z$bi),
                anindex=as.double(z$anindex),
                andirection=as.double(z$andirection),
                det=as.double(z$det),
                as.double(Bcov),
                as.double(solvebtb),
                as.double(projectmat),
                as.double(sigma2),
                as.double(z$sigma2hat),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.integer(ngrad),
                as.double(hakt),
                as.double(zext),
                as.double(rho),
                as.double(lambda0),
                theta=double(6*n),
                s0hat=double(n),
                sigma2hat=double(n),
                double(ngrad),
                as.double(eps),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","bi","anindex","andirection","det","s0hat","sigma2hat")]
     if(hakt<hsig){
        eta <- (hsig^3 - hakt^3)/hsig^3
        z$sigma2hat <- eta*sigma2hat+(1-eta)*z$sigma2hat
     }
     dim(z$s0hat) <- dim(z$bi) <- dim(z$anindex) <- dim(z$det) <- dim(z$sigma2hat) <- dimy[-1]
     dim(z$theta) <- dimy
     dim(z$andirection) <- c(3,dimy[-1]) 
     if(graph){
     class(z) <- "dti"
     img<-z$theta[1,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dxx: mean",signif(mean(img),3),"max",signif(max(img),3)))
     img<-z$theta[2,,,slice]
     show.image(make.image(img))
     title(paste("Dxy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[3,,,slice]
     show.image(make.image(img))
     title(paste("Dxz: min",signif(min(img),3),"max",signif(max(img),3)))
     show.image(make.image(z$anindex[,,slice]))
     title(paste("Anisotropy index  range:",signif(min(z$anindex),3),"-",
                  signif(max(z$anindex),3)))
     img<-z$theta[4,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: mean",signif(mean(img),3),"max",signif(max(img),3)))
     img<-z$theta[5,,,slice]
     show.image(make.image(img))
     title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
#     andir2.image(z,slice,quant=quant,minanindex=minanindex)
     plot(c(0,1),c(0,1))
     title(paste("Directions (h=",signif(hakt,3),"), slice",slice))
     ni<-z$bi[,,slice]*z$sigma2hat[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights  mean=",signif(mean(z$bi*z$sigma2hat),3)))
     img<-z$theta[6,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: mean",signif(mean(img),3),"max",signif(max(img),3)))
     }
     cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex,c(.5, .75, .9, .95, 1)),3),"\n")
     hakt0<-hakt
     hakt <- hakt*hincr
  }
  dim(z$s0hat) <- c(n1,n2,n3)
#  z <- list(theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,
#            ddim0=ddim0,xind=xind,yind=yind,zind=zind,InvCov=Bcov,s0hat=z$s0hat,call=args)

  invisible(new("dtiTensor",
                list(theta = z$theta, sigma = z$sigma2hat, scorr = scorr),
                btb   = btb,
                ngrad = ngrad, # = dim(btb)[2]
                ddim  = as.integer(ddim),
                ddim0 = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                source= source)
            )

})

setMethod("dti.smooth", "dtiTensor", function(object,method="riemann",hmax=5,lambda=20,rho=1,graph=FALSE,
                                            slice=NULL,quant=.8,minanindex=NULL,zext=1){
  if (method == "riemann") {
    rdtianiso(object,hmax=hmax,lambda=lambda,rho=rho,graph=graph,
              slice=slice,quant=quant,minanindex=minanindex,zext=zext)
  } else {
    dtianiso(object,hmax=hmax,lambda=lambda,rho=rho,graph=graph,
              slice=slice,quant=quant,minanindex=minanindex,zext=zext)
  }
})

dtianiso <- function(dtobject,hmax=5,lambda=20,rho=1,graph=FALSE,slice=NULL,quant=.8,minanindex=NULL,zext=1){
  cat("Smoothing DT with Euclidean Metric\n")
  if (graph) {
    adimpro <- require(adimpro)
    if (!adimpro) cat("No graphical output! Install package adimpro from CRAN!\n")
    graph <- graph & adimpro
  }
  args <- match.call()
  btb<-dtobject@btb
  Bcov <- btb%*%t(btb)
  y <- dtobject$theta
  sigma2 <- dtobject$sigma
  scorr <- dtobject$scorr
  h0 <- dtobject$bw
  mask <- array(1,dim=dtobject@ddim)
  ddim <- dtobject@ddim
  ddim0 <- dtobject@ddim0
  xind <- dtobject@xind
  yind <- dtobject@yind
  zind <- dtobject@zind
  source <- dtobject@source
  rm(dtobject)
  gc()
  dimy <- dim(y)
  if(length(dimy)!=4||dimy[1]!=6) stop("y does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  sigma2[sigma2<=mean(sigma2)*1e-5]<- mean(sigma2)*1e-5
  z <- .Fortran("projdt",
                as.double(y),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                theta=double(6*n),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                mask=as.logical(mask),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","anindex","andirection","det","mask")]
  y <- array(z$theta,dimy)
  z$bi <- 1/sigma2
  dim(z$theta) <- dimy
  dim(z$anindex) <-dim(z$det) <-dim(z$mask) <- dimy[-1]
  z$mask <- array(z$mask,dimy[-1])&mask
  dim(z$andirection) <- c(3,dimy[-1]) 
#
#  initial state for h=1
#
  if(graph){
    oldpar <- par(mfrow=c(3,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
    on.exit(par(oldpar))
    if(is.null(slice)) slice<-n3%/%2
    class(z) <- "dti"
    img<-z$theta[1,,,slice]
    show.image(make.image(65535*img/max(img)))
    title(paste("Dxx: min",signif(min(img),3),"max",signif(max(img),3)))
    img<-z$theta[2,,,slice]
    show.image(make.image(img))
    title(paste("Dxy: min",signif(min(img),3),"max",signif(max(img),3)))
    img<-z$theta[3,,,slice]
    show.image(make.image(img))
    title(paste("Dxz: min",signif(min(img),3),"max",signif(max(img),3)))
    show.image(make.image(z$anindex[,,slice]))
    title(paste("Anisotropy index  range:",signif(min(z$anindex[z$mask]),3),"-",
                signif(max(z$anindex[z$mask]),3)))
    img<-z$theta[4,,,slice]
    show.image(make.image(65535*img/max(img)))
    title(paste("Dyy: min",signif(min(img),3),"max",signif(max(img),3)))
    img<-z$theta[5,,,slice]
    show.image(make.image(img))
    title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
    andir.image(z,slice,quant=quant,minanindex=minanindex)
    title(paste("Directions (h=1), slice",slice))
    ni<-z$bi[,,slice]*sigma2[,,slice]
    show.image(make.image(65535*ni/max(ni)))
    title(paste("sum of weights  mean=",signif(mean(z$bi[z$mask]*sigma2[z$mask]),3)))
    img<-z$theta[6,,,slice]
    show.image(make.image(65535*img/max(img)))
    title(paste("Dzz: min",signif(min(img),3),"max",signif(max(img),3)))
  }
  hincr <- 1.25^(1/3)
  hakt0 <- 1
  hakt <- hincr
  lambda0 <- lambda
  while( hakt <= hmax) {
    if (any(h0 >= 0.25)) {
      corrfactor <- Spatialvar.gauss(hakt0/0.42445/4,h0,3) /
        Spatialvar.gauss(h0,1e-5,3) /
          Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
      lambda0 <- lambda * corrfactor
      cat("Correction factor for spatial correlation",signif(corrfactor,3),"\n")
    }
    z <- .Fortran("awsdti",
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
    dim(z$bi) <- dim(z$anindex) <- dim(z$det) <- dim(z$mask) <- dimy[-1]
    dim(z$theta) <- dimy
    dim(z$andirection) <- c(3,dimy[-1]) 
    if(graph){
      class(z) <- "dti"
      img<-z$theta[1,,,slice]
      show.image(make.image(65535*img/max(img)))
      title(paste("Dxx: min",signif(min(img),3),"max",signif(max(img),3)))
      img<-z$theta[2,,,slice]
      show.image(make.image(img))
      title(paste("Dxy: min",signif(min(img),3),"max",signif(max(img),3)))
      img<-z$theta[3,,,slice]
      show.image(make.image(img))
      title(paste("Dxz: min",signif(min(img),3),"max",signif(max(img),3)))
      show.image(make.image(z$anindex[,,slice]))
      title(paste("Anisotropy index  range:",signif(min(z$anindex[z$mask]),3),"-",
                  signif(max(z$anindex[z$mask]),3)))
      img<-z$theta[4,,,slice]
      show.image(make.image(65535*img/max(img)))
      title(paste("Dyy: min",signif(min(img),3),"max",signif(max(img),3)))
      img<-z$theta[5,,,slice]
      show.image(make.image(img))
      title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
      andir.image(z,slice,quant=quant,minanindex=minanindex)
      title(paste("Directions (h=",signif(hakt,3),"), slice",slice))
      ni<-z$bi[,,slice]*sigma2[,,slice]
      show.image(make.image(65535*ni/max(ni)))
      title(paste("sum of weights  mean=",signif(mean(z$bi[z$mask]*sigma2[z$mask]),3)))
      img<-z$theta[6,,,slice]
      show.image(make.image(65535*img/max(img)))
      title(paste("Dyy: min",signif(min(img),3),"max",signif(max(img),3)))
    }
    cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex[z$mask],c(.5, .75, .9, .95, 1)),3),"\n")
    hakt0<-hakt
    hakt <- hakt*hincr
  }
#  z <- list(theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,mask=z$mask,InvCov=Bcov,call=args)
#  class(z) <- "dti"

  invisible(new("dtiTensor",
                list(theta = z$theta, sigma = z$sigma2hat, scorr = scorr),
                btb   = btb,
                ngrad = ngrad, # = dim(btb)[2]
                ddim  = as.integer(ddim),
                ddim0 = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                source= source)
            )
}


rdtianiso <- function(dtobject,hmax=5,lambda=20,rho=1,graph=FALSE,slice=NULL,quant=.8,minanindex=NULL,zext=1){
  cat("Smoothing DT with Riemann Tensor Metric\n")
  if (graph) {
    adimpro <- require(adimpro)
    if (!adimpro) cat("No graphical output! Install package adimpro from CRAN!\n")
    graph <- graph & adimpro
  }
  args <- match.call()
  btb<-dtobject@btb
  Bcov <- btb%*%t(btb)
  y <- dtobject$theta
  sigma2 <- dtobject$sigma
  scorr <- dtobject$scorr
  h0 <- dtobject$bw
  mask <- array(1,dim=dtobject@ddim)
  ddim <- dtobject@ddim
  ddim0 <- dtobject@ddim0
  xind <- dtobject@xind
  yind <- dtobject@yind
  zind <- dtobject@zind
  source <- dtobject@source
  rm(dtobject)
  gc()
  dimy <- dim(y)
  if(length(dimy)!=4||dimy[1]!=6) stop("y does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  sigma2[sigma2<=mean(sigma2)*1e-5]<- mean(sigma2)*1e-5
  z <- .Fortran("projdt",
                as.double(y),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                theta=double(6*n),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                mask=as.logical(mask),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","anindex","andirection","det","mask")]
  y <- array(z$theta,dimy)
  z$bi <- 1/sigma2
  dim(z$theta) <- dimy
  dim(z$anindex) <-dim(z$det) <-dim(z$mask) <- dimy[-1]
  z$mask <- array(z$mask,dimy[-1])&mask
  dim(z$andirection) <- c(3,dimy[-1]) 
#
#  initial state for h=1
#
  if(graph){
     oldpar <- par(mfrow=c(3,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
     on.exit(par(oldpar))
     if(is.null(slice)) slice<-n3%/%2
     class(z) <- "dti"
     img<-z$theta[1,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dxx: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[2,,,slice]
     show.image(make.image(img))
     title(paste("Dxy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[3,,,slice]
     show.image(make.image(img))
     title(paste("Dxz: min",signif(min(img),3),"max",signif(max(img),3)))
     show.image(make.image(z$anindex[,,slice]))
     title(paste("Anisotropy index  range:",signif(min(z$anindex[z$mask]),3),"-",
                  signif(max(z$anindex[z$mask]),3)))
     img<-z$theta[4,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[5,,,slice]
     show.image(make.image(img))
     title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
     andir.image(z,slice,quant=quant,minanindex=minanindex)
     title(paste("Directions (h=1), slice",slice))
     ni<-z$bi[,,slice]*sigma2[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights  mean=",signif(mean(z$bi[z$mask]*sigma2[z$mask]),3)))
     img<-z$theta[6,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dzz: min",signif(min(img),3),"max",signif(max(img),3)))
  }
  hincr <- 1.25^(1/3)
  hakt0 <- 1
  hakt <- hincr
  lambda0 <- lambda
  while( hakt <= hmax) {
    if (any(h0 >= 0.25)) {
       corrfactor <- Spatialvar.gauss(hakt0/0.42445/4,h0,3) /
        Spatialvar.gauss(h0,1e-5,3) /
        Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
        lambda0 <- lambda * corrfactor
     cat("Correction factor for spatial correlation",signif(corrfactor,3),"\n")
}
     usize<-(2*ceiling(hakt+1))^3/zext
     mask[z$anindex<0.1]<-FALSE
     z <- .Fortran("rawsdti",
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
                double(usize),# array for weights
                double(9*usize),# array for complete diffusion tensors with pos. weights
                as.integer(usize),# upper bound for size of U(x)
                double(9*usize),# array for Remannian Log Maps 
                DUP=FALSE,
                PACKAGE="dti")[c("theta","bi","anindex","andirection","det","mask")]
     dim(z$bi) <- dim(z$anindex) <- dim(z$det) <- dim(z$mask) <- dimy[-1]
     dim(z$theta) <- dimy
     dim(z$andirection) <- c(3,dimy[-1]) 
     if(graph){
     
     class(z) <- "dti"
     img<-z$theta[1,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dxx: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[2,,,slice]
     show.image(make.image(img))
     title(paste("Dxy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[3,,,slice]
     show.image(make.image(img))
     title(paste("Dxz: min",signif(min(img),3),"max",signif(max(img),3)))
     show.image(make.image(z$anindex[,,slice]))
     title(paste("Anisotropy index  range:",signif(min(z$anindex[z$mask]),3),"-",
                  signif(max(z$anindex[z$mask]),3)))
     img<-z$theta[4,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[5,,,slice]
     show.image(make.image(img))
     title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
     andir.image(z,slice,quant=quant,minanindex=minanindex)
     title(paste("Directions (h=",signif(hakt,3),"), slice",slice))
     ni<-z$bi[,,slice]*sigma2[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights  mean=",signif(mean(z$bi[z$mask]*sigma2[z$mask]),3)))
     img<-z$theta[6,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: min",signif(min(img),3),"max",signif(max(img),3)))
     }
     cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex[z$mask],c(.5, .75, .9, .95, 1)),3),"\n")
     hakt0<-hakt
     hakt <- hakt*hincr
  }
#z <- list(theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,mask=z$mask, ddim0=ddim0,xind=xind,yind=yind,zind=zind,InvCov=Bcov,call=args)
#class(z) <- "dti"

  invisible(new("dtiTensor",
                list(theta = z$theta, sigma = NULL, scorr = scorr),
                btb   = btb,
                ngrad = ngrad, # = dim(btb)[2]
                ddim  = as.integer(ddim),
                ddim0 = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                source= source)
            )
}

