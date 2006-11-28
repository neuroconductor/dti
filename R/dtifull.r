create2.dti <- function(gradient,imagefile,ddim,xind=NULL,yind=NULL,zind=NULL){
if(dim(gradient)[2]==3)  gradient<-t(gradient)
if(dim(gradient)[1]!=3)  stop("Not a valid gradient matrix")
ngrad <- dim(gradient)[2]
if(!(file.exists(imagefile))) stop("Image file does not exist")
zz<-file(imagefile,"rb")
s0 <- readBin(zz,"integer",prod(ddim),2,FALSE)
si <- readBin(zz,"integer",prod(ddim)*ngrad,2,FALSE)
close(zz)
cat("Data successfully read \n")
if(is.null(xind)) xind<-1:ddim[1]
if(is.null(yind)) yind<-1:ddim[2]
if(is.null(zind)) zind<-1:ddim[3]
dim(s0) <- ddim
s0 <- s0[xind,yind,zind]
dim(si) <- c(ddim,ngrad)
si <- si[xind,yind,zind,]
ddim0 <- ddim
ddim <- dim(s0)
dim(s0)<-dim(si)<-NULL
ttt <- -log(si/s0)
ttt[is.na(ttt)] <- 0
ttt[(ttt==Inf)] <- 0
ttt[(ttt==-Inf)] <- 0
n <- prod(ddim)
dim(ttt) <- c(n,ngrad)
ttt<-t(ttt)
cat("Data transformation completed \n")
btb <- matrix(0,6,ngrad)
btb[1,]<-gradient[1,]*gradient[1,]
btb[4,]<-gradient[2,]*gradient[2,]
btb[6,]<-gradient[3,]*gradient[3,]
btb[2,]<-2*gradient[1,]*gradient[2,]
btb[3,]<-2*gradient[1,]*gradient[3,]
btb[5,]<-2*gradient[2,]*gradient[3,]
btbsvd <- svd(btb)
solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
theta <- solvebtb%*% ttt
cat("Diffusion tensors generated \n")
res <- ttt - t(btb) %*% theta
mres2 <- res[1,]^2
for(i in 2:ngrad) mres2 <- mres2 + res[i,]^2
sigma2 <- array(mres2/(ngrad-6),ddim)
cat("Variance estimates generated \n")
rm(mres2)
gc()
z<-list(theta=array(theta,c(6,ddim)),sigma2=sigma2,btb=btb,
        solvebtb=solvebtb,scorr=c(0,0),s0=array(s0,ddim),
        ddim=ddim,ddim0=ddim0,xind=xind,yind=yind,zind=zind,
        ngrad=ngrad,file=imagefile,si=si,res=res)
class(z) <- "dti"
invisible(z)
}

getscorr2 <- function(dtobject,level=0){
if(!("dti" %in% class(dtobject))) stop("Not an dti-object")
ddim <- dtobject$ddim
mask <- dtobject$s0>level
n <- prod(ddim)
ngrad <- dtobject$ngrad
scorr <- c(0,0)
res <- dtobject$res
dim(res) <- c(ngrad,ddim)
res <- aperm(res,c(2:4,1))
dim(res) <- c(n,ngrad)
res1 <- as.vector(res[as.vector(mask),])
scorr[1] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
cat("correlation in x-direction",signif(scorr[1],3),"\n")
dim(res) <- c(ddim,ngrad)
res <- aperm(res,c(2,1,3,4))
dim(res) <- c(n,ngrad)
res1 <- as.vector(res[as.vector(aperm(mask,c(2,1,3))),])
scorr[2] <- mean(res1[-1]*res1[-length(res1)])/var(res1)
cat("correlation in y-direction",signif(scorr[2],3),"\n")
dtobject$scorr <- scorr
dtobject$res <- NULL
invisible(dtobject)
}






dtianiso2<-function(dtobject,hmax=5,lambda=20,rho=1,graph=FALSE,slice=NULL,quant=.8,minanindex=NULL,zext=1){
if(!("dti" %in% class(dtobject))) stop("Not an dti-object")
  args <- match.call()
  btb<-dtobject$btb
  Bcov <- btb%*%t(btb)
  y <- dtobject$theta
  sigma2 <- dtobject$sigma2
  scorr <- dtobject$scorr
  si <- dtobject$si
  s0 <- dtobject$s0
  ngrad <- dtobject$ngrad
  solvebtb <- dtobject$solvebtb
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
  z <- .Fortran("projdt2",
                as.double(y),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                theta=double(6*n),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","anindex","andirection","det")]
  y <- array(z$theta,dimy)
  z$bi <- 1/sigma2
  dim(z$theta) <- dimy
  dim(z$anindex) <-dim(z$det) <- dimy[-1]
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
     title(paste("Anisotropy index  range:",signif(min(z$anindex),3),"-",
                  signif(max(z$anindex),3)))
     img<-z$theta[4,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[5,,,slice]
     show.image(make.image(img))
     title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
     andir2.image(z,slice,quant=quant,minanindex=minanindex)
     title(paste("Directions (h=1), slice",slice))
     ni<-z$bi[,,slice]*sigma2[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights  mean=",signif(mean(z$bi*sigma2),3)))
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
                as.double(sigma2),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.integer(ngrad),
                as.double(hakt),
                as.double(zext),
                as.double(rho),
                as.double(lambda0),
                theta=double(6*n),
                double(ngrad),
                DUP=FALSE,
                PACKAGE="dti")[c("theta","bi","anindex","andirection","det")]
     dim(z$bi) <- dim(z$anindex) <- dim(z$det) <- dimy[-1]
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
     title(paste("Anisotropy index  range:",signif(min(z$anindex),3),"-",
                  signif(max(z$anindex),3)))
     img<-z$theta[4,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: min",signif(min(img),3),"max",signif(max(img),3)))
     img<-z$theta[5,,,slice]
     show.image(make.image(img))
     title(paste("Dyz: min",signif(min(img),3),"max",signif(max(img),3)))
     andir2.image(z,slice,quant=quant,minanindex=minanindex)
     title(paste("Directions (h=",signif(hakt,3),"), slice",slice))
     ni<-z$bi[,,slice]*sigma2[,,slice]
     show.image(make.image(65535*ni/max(ni)))
     title(paste("sum of weights  mean=",signif(mean(z$bi*sigma2),3)))
     img<-z$theta[6,,,slice]
     show.image(make.image(65535*img/max(img)))
     title(paste("Dyy: min",signif(min(img),3),"max",signif(max(img),3)))
     }
     cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex,c(.5, .75, .9, .95, 1)),3),"\n")
     hakt0<-hakt
     hakt <- hakt*hincr
  }
z <- list(theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,
          ddim0=ddim0,xind=xind,yind=yind,zind=zind,InvCov=Bcov,call=args)
class(z) <- "dti"
invisible(z)
}
andir2.image <- function(dtobject,slice=1,method=1,quant=0,minanindex=NULL,show=TRUE,...){
if(!("dti" %in% class(dtobject))) stop("Not an dti-object")
if(is.null(dtobject$anindex)) stop("No anisotropy index yet")
adimpro <- require(adimpro)
anindex <- dtobject$anindex
dimg <- dim(anindex)[1:2]
if(is.null(slice)) slice <- 1
anindex <- anindex[,,slice]
andirection <- dtobject$andirection[,,,slice]
anindex[anindex>1]<-0
anindex[anindex<0]<-0
dim(andirection)<-c(3,prod(dimg))
if(is.null(minanindex)) minanindex <- quantile(anindex,quant)
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
andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
dim(andirection)<-c(dimg,3)
if(adimpro) {
andirection <- make.image(andirection)
if(show) show.image(andirection,...)
} else if(show) {
dim(anindex) <- dimg
image(anindex,...)
}
invisible(andirection)
} 
