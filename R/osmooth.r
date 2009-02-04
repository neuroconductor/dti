risk <- function(andir1,dtobj2,minvalue){
andir2 <- extract(dtobj2,"andir")$andir
fa <- extract(dtobj2,"fa")$fa
mask <- extract(dtobj2,"mask")$mask
mask <- mask&fa>minvalue
n <- prod(dtobj2@ddim)
dim(andir1) <- dim(andir2) <- c(3,n)
sqrt(1-mean((c(1,1,1)%*%(andir1[,mask]*andir2[,mask]))^2))
}
#
#   Nonlinear regression; regularized version according to Koay et. al. (2006)
#   this is also based on an a statistical penalty defined using log-likelihood difference
#
dtireg.osmooth <- function(object,hmax=5,hinit=1,lambda=30,tau=10,graph=FALSE,slice=NULL,quant=.8,
                         minanindex=NULL,hsig=2.5,lseq=NULL,varmethod="residuals",rician=FALSE,niter=5,varmodel="local",result="Tensor"){
#
#     lambda and lseq adjusted for alpha=0.2
#
  if(rician) {
     warning("Rician bias correction not yet implemented for this method")
     rician <- FALSE
  }
  if(!is.null(object$ni)){
     warning("DWI object has been smoothed already, smoothing omitted")
     return(if(result=="Tensor") dtiTensor(object,method="nonlinear",varmethod=varmethod,varmodel=varmodel) else object)
  }
  eps <- 1e-6
  maxnw <- 10000
  if (graph) {
    adimpro <- require(adimpro)
    if (!adimpro) cat("No graphical output! Install package adimpro from CRAN!\n")
    graph <- graph & adimpro
  }
  args <- sys.call(-3)
  args <- c(object@call,args)
  dtobject <- dtiTensor(object,method="nonlinear",varmethod=varmethod,varmodel=varmodel)
  scale <- dtobject@scale
  mask <- dtobject@mask
  th0 <- dtobject@th0
  D <- dtobject@D
  scorr <- dtobject@scorr
  h0 <- dtobject@bw
  rm(dtobject)
  gc()
  cat("Corresponding bandwiths for specified correlation:",h0,"\n")
  s0ind <- object@s0ind
  ngrad <- object@ngrad
  ddim0 <- object@ddim0
  ddim <- object@ddim
  sdcoef <- object@sdcoef
  z <- .Fortran("outlier",
                as.integer(object@si),
                as.integer(prod(ddim)),
                as.integer(ngrad),
                as.logical((1:ngrad)%in%s0ind),
                as.integer(length(s0ind)),
                si=integer(prod(ddim)*ngrad),
                index=integer(prod(ddim)),
                lindex=integer(1),
                DUPL=FALSE,
                PACKAGE="dti")[c("si","index","lindex")]
  si <- array(z$si,c(ddim,ngrad))
  index <- if(z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  gc()
  si <- aperm(si,c(4,1:3))
  xind <- object@xind
  yind <- object@yind
  zind <- object@zind
  source <- object@source
  btb <- object@btb
  voxelext <- object@voxelext
  if(is.null(voxelext)) vext <- c(1,1,1) else vext <- voxelext/min(voxelext)
  gc()
  dimy <- dim(D)
  if(length(dimy)!=4||dimy[1]!=6) stop("D does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  z <- .Fortran("projdt2",
                as.double(D),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                D=double(6*n),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                as.double(eps),
                DUP=FALSE,
                PACKAGE="dti")[c("D","anindex","andirection","det")]
  dim(z$D) <- dimy
  z$th0 <- th0
  dim(z$anindex) <-dim(z$det) <- dimy[-1]
  dim(z$andirection) <- c(3,dimy[-1]) 
  z$sihat <- si
  low <- sdcoef[1]+sdcoef[3]*sdcoef[2]
  high <- sdcoef[1]+sdcoef[4]*sdcoef[2]
  z$varinv <- 1/pmin(high,pmax(low,sdcoef[1]+si*sdcoef[2]))^2
#
#  initial state for h=1
#
  if(graph){
     oldpar <- par(mfrow=c(3,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
     on.exit(par(oldpar))
     if(is.null(slice)) slice<-n3%/%2
     class(z) <- "dti"
     img<-z$D[1,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dxx: mean",signif(mean(z$D[1,,,][mask]),3),"max",signif(max(z$D[1,,,][mask]),3)))
     img<-z$D[2,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxy: min",signif(min(z$D[2,,,][mask]),3),"max",signif(max(z$D[2,,,][mask]),3)))
     img<-z$D[3,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxz: min",signif(min(z$D[3,,,][mask]),3),"max",signif(max(z$D[3,,,][mask]),3)))
     show.image(make.image(z$anindex[,,slice]),xaxt="n",yaxt="n")
     title(paste("Anisotropy index  range:",signif(min(z$anindex[mask]),3),"-",
                  signif(max(z$anindex[mask]),3)))
     img<-z$D[4,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dyy: mean",signif(mean(z$D[4,,,][mask]),3),"max",signif(max(z$D[4,,,][mask]),3)))
     img<-z$D[5,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dyz: min",signif(min(z$D[5,,,][mask]),3),"max",signif(max(z$D[5,,,][mask]),3)))
     andir2.image(z,slice,quant=quant,minanindex=minanindex,xaxt="n",yaxt="n")
     title(paste("Directions (h=1), slice",slice))
     ni <- array(1,dimy[-1])*as.integer(mask)
     show.image(make.image((65535*ni/max(ni))[,,slice]),xaxt="n",yaxt="n")
     title(paste("sum of weights  mean=",signif(1,3)))
     img<-z$D[6,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dzz: mean",signif(mean(z$D[6,,,][mask]),3),"max",signif(max(z$D[6,,,][mask]),3)))
  }
  hincr <- 1.25^(1/3)
  maxvol <- getvofh(hmax,c(1,0,0,1,0,1),vext)
  kstar <- as.integer(log(maxvol)/log(1.25)) 
  if(is.null(hinit)){
  hakt0 <- 1
  hakt <- hincr
  hinit <- 1
  } else {
  hakt0 <- max(1,hinit/hincr)
  hakt <- hinit
  }
  steps <- kstar+1

  # define lseq
  if (is.null(lseq)) {
# this is optimized for lkern="Gaussian" such that alpha approx 0.04 -- 0.1 and probability of separated points is approx. 1e-4
    lseq <- c(1.5,.9,.8,.85,.8,0.85,.85,0.95,1.25,1.15)# alpha=0.2    abs. deviation in S0
  }
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
  k <- 1
  lambda0 <- lambda
  while(k <= kstar) {
      hakt0 <- gethani(1,10,1.25^(k-1),c(1,0,0,1,0,1),vext,1e-4)
      hakt <- gethani(1,10,1.25^k,c(1,0,0,1,0,1),vext,1e-4)
    if (any(h0 >= 0.25)) {
       corrfactor <- Spatialvar.gauss(hakt0/0.42445/4,h0,3) /
       Spatialvar.gauss(h0,1e-5,3) /
       Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
       lambda0 <- lambda0 * corrfactor
       cat("Correction factor for spatial correlation",signif(corrfactor,3),"\n")
    }
    z <- .Fortran("osmdti",
                    as.integer(si),
                    sihat=as.integer(z$sihat), # needed for statistical penalty
                    varinv=as.double(z$varinv), # needed for statistical penalty
                    double(ngrad*n),# array for predicted Si's from the tensor model 
                    as.integer(ngrad),
                    as.integer(n1),
                    as.integer(n2),
                    as.integer(n3),
                    as.logical(mask),
                    as.double(btb),
                    as.double(sdcoef),
                    as.double(z$th0),
                    th0=double(n),
                    as.double(z$D),
                    D=double(6*n),
                    anindex=as.double(z$anindex),
                    andirection=as.double(z$andirection),
                    det=as.double(z$det),
                    as.double(1.25^k),
                    as.integer(niter),
                    as.double(vext),
                    as.double(lambda0),
                    as.double(tau),
                    double(ngrad),#wijk
                    double(ngrad),#swsi
                    double(ngrad),#swsi2
                    double(ngrad),#sw
                    double(ngrad),#sw2
                    double(ngrad),#oibd
                    double(ngrad),#F
                    double(ngrad),#var
                    double(ngrad),#sig2bi
                    as.logical((1:ngrad)%in%s0ind),
                    as.double(eps),
                    DUP=FALSE,PACKAGE="dti")[c("th0","D","anindex","andirection","det","sihat","varinv")] 
     dim(z$th0) <- dim(z$anindex) <- dim(z$det) <- dimy[-1]
     dim(z$D) <- dimy
     dim(z$andirection) <- c(3,dimy[-1]) 
     if(any(is.na(z$th0))){
        indna <- is.na(z$th0)
        cat("found ",sum(indna),"NA's\n")
        z$th0[indna] <- th0[indna]
        cat("th0",th0[indna])
        dim(z$D) <- dim(D) <- c(6,n1*n2*n3)
        z$D[,indna] <- D[,indna]
        dim(z$D) <- dim(D) <- c(6,n1,n2,n3)
        mask[indna] <- FALSE
     }
##
##  set mask to FALSE on voxel where we observe extreme 
##  values of z$det
##  this is an indication for numerical problems caused
##  e.g. by registration artifacts
##
     det95 <- quantile(z$det[mask],.95)
     if(any(z$det[mask]>1e3*det95)){
        n <- n1*n2*n3
        indna <- (1:n)[z$det>1e3*det95]
        cat("found ",length(indna),"voxel with extreme derterminat 's, keep initial estimates and do not use these voxel\n")
        z$th0[indna] <- th0[indna]
        dim(z$D) <- dim(D) <- c(6,n1*n2*n3)
        z$D[,indna] <- D[,indna]
        dim(z$D) <- dim(D) <- c(6,n1,n2,n3)
        dim(z$sihat)  <- dim(si) <- c(ngrad,n1*n2*n3)
        z$sihat[,indna] <- si[,indna]
        dim(z$sihat) <- dim(si) <- c(ngrad,n1,n2,n3)
        z$det[indna] <- 0
        mask[indna] <- FALSE
     }
     if(graph){
     class(z) <- "dti"
     img<-z$D[1,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dxx: mean",signif(mean(z$D[1,,,][mask]),3),"max",signif(max(z$D[1,,,][mask]),3)))
     img<-z$D[2,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxy: min",signif(min(z$D[2,,,][mask]),3),"max",signif(max(z$D[2,,,][mask]),3)))
     img<-z$D[3,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxz: min",signif(min(z$D[3,,,][mask]),3),"max",signif(max(z$D[3,,,][mask]),3)))
     show.image(make.image(z$anindex[,,slice]),xaxt="n",yaxt="n")
     title(paste("Anisotropy index  range:",signif(min(z$anindex[mask]),3),"-",
                  signif(max(z$anindex[mask]),3)))
     img<-z$D[4,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dyy: mean",signif(mean(z$D[4,,,][mask]),3),"max",signif(max(z$D[4,,,][mask]),3)))
     img<-z$D[5,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dyz: min",signif(min(z$D[5,,,][mask]),3),"max",signif(max(z$D[5,,,][mask]),3)))
     andir2.image(z,slice,quant=quant,minanindex=minanindex,xaxt="n",yaxt="n")
     title(paste("Directions (h=",signif(hakt,3),"), slice",slice))
     dim(z$varinv) <- c(ngrad,n1,n2,n3)
     ni <- apply(z$varinv[,,,slice],2:3,mean)*mask[,,slice]
     show.image(make.image(65535*ni/max(ni)),xaxt="n",yaxt="n")
     title(paste("sum of weights  mean=",signif(mean(z$varinv[mask]),3)))
     img<-z$D[6,,,slice]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dzz: mean",signif(mean(z$D[6,,,][mask]),3),"max",signif(max(z$D[6,,,][mask]),3)))
     }
     cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex[mask],c(.5, .75, .9, .95, 1)),3),"\n")
     cat("risk=",risk(z$andirection,dt0,0.55),risk(z$andirection,dt0,0.65),risk(z$andirection,dt0,0.75),risk(z$andirection,dt0,0.85),"\n")
    k <- k+1
     lambda0 <- lambda 
     gc()
  }
  dimsi <- dim(si)
  rm(si)
  gc()
  if(result=="Tensor"){
  cat("prepare final dtiTensor object",date(),"\n")
  } else   cat("prepare final smoothed dtiData object",date(),"\n")
  if(result=="Tensor") invisible(new("dtiTensor",
                list(ni=z$varinv),
                call = args,
                D = z$D,
                th0 = z$th0,
                sigma = array(0,c(1,1,1)),
                scorr = scorr,
                bw = h0,
                mask = mask,
                btb   = btb,
                hmax  = hmax,
                ngrad = ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = as.integer(ddim),
                ddim0 = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                voxelext = object@voxelext,
                source= object@source,
                outlier = index,
                scale = scale,
                method= "nonlinear")
            ) else invisible(new("dtiData",
                list(ni=z$varinv),
                call = args,
                si = aperm(array(as.integer(z$sihat),dimsi),c(2:4,1)),
                sdcoef = sdcoef,
                btb    = btb,
                ngrad  = ngrad, # = dim(btb)[2]
                s0ind  = object@s0ind, # indices of s0 images
                replind = object@replind, # replications in gradient design
                ddim   = as.integer(ddim),
                ddim0  = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                voxelext = object@voxelext,
                level  = object@level,
                orientation = object@orientation ,
                source= object@source)
            )
}


