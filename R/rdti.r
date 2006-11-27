rdtianiso<-function(dtobject,hmax=5,lambda=20,rho=1,graph=FALSE,slice=NULL,quant=.8,minanindex=NULL,zext=1){
if(!("dti" %in% class(dtobject))) stop("Not an dti-object")
  args <- match.call()
  btb<-dtobject$btb
  Bcov <- btb%*%t(btb)
  y <- dtobject$theta
  sigma2 <- dtobject$sigma2
  scorr <- dtobject$scorr
  mask <- dtobject$mask
  ddim0 <- dtobject$ddim0
  xind <- dtobject$xind
  yind <- dtobject$yind
  zind <- dtobject$zind
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
     require(adimpro)
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
  if (max(scorr)>0) {
    h0 <- numeric(length(scorr))
    for (i in 1:length(h0)) h0[i] <- get.bw.gauss(scorr[i],interv=2)
    if (length(h0)<2) h0 <- rep(h0[1],2)
    h0 <- c(h0,1e-5)
# no spatial correlation iz z-direction
    cat("Corresponding bandwiths for specified correlation:",h0,"\n")
  }
  hincr <- 1.25^(1/3)
  hakt0 <- 1
  hakt <- hincr
  lambda0 <- lambda
  while( hakt <= hmax) {
    if (scorr[1]>=0.1) {
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
z <- list(theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,mask=z$mask,
          ddim0=ddim0,xind=xind,yind=yind,zind=zind,InvCov=Bcov,call=args)
class(z) <- "dti"
invisible(z)
}

test <- function(mat1,mat2){
rlm<-.Fortran("rlogmap",
             as.double(mat1),
             as.double(mat2),
             as.integer(1),
             rlm=double(9),
             DUP=FALSE,
             PACKAGE="dti")$rlm
rem <- .Fortran("rexpmap",
             as.double(mat1),
             as.double(rlm),
             rem=double(9),
             integer(1),
             DUP=FALSE,
             PACKAGE="dti")$rem
print(mat1)
print(mat2)
print(matrix(rem,3,3))
}