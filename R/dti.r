#
#
#

setMethod("print", "dti",
function(x){
    cat("  DTI object of class", class(x),"\n")
    cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(x@ngrad, collapse="x"), "\n")
    cat("  Filename             :", x@source, "\n")
    cat("  Slots                :\n")
    print(slotNames(x))
    invisible(NULL)
})
setMethod("summary", "dti",
function(object){
    cat("  DTI object of class", class(object),"\n")
    cat("  Filename             :", object@source, "\n")
    cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
    cat("  Number of Gradients  :", paste(object@ngrad, collapse="x"), "\n")
    cat("  Voxel extensions     :", paste(object@voxelext, collapse="x"), "\n")
    cat("  Ind. S0-Images       :", paste(object@s0ind, collapse="x"), "\n")
    if(class(object)=="dtiTensor"){
    cat("  Voxel in mask        :", paste(sum(object@mask), collapse="x"), "\n")
    cat("  Spatial smoothness   :", paste(signif(object@bw,3), collapse="x"), "\n")
    cat("  mean variance        :", paste(signif(mean(object@sigma[object@mask]),3), collapse="x"), "\n")
    cat("  hmax                 :", paste(object@hmax, collapse="x"), "\n")
}
    cat("\n")
    invisible(NULL)
})

setMethod("plot", "dtiTensor", function(x, y, slice=1, view="axial", quant=0, minanindex=NULL, contrast.enh=1,what="FA", qrange=c(.01,.99),xind=NULL,yind=NULL,zind=NULL, mar=c(2,2,2,.2),mgp=c(2,1,0),...) {
  if(is.null(x@D)) cat("No diffusion tensor yet")
  adimpro <- require(adimpro)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  if (view == "sagittal") {
    D <- x@D[,slice,yind,zind]
    mask <- x@mask[slice,yind,zind]
    n1 <- length(yind)
    n2 <- length(zind)
  } else if (view == "coronal") {
    D <- x@D[,xind,slice,zind]
    mask <- x@mask[xind,slice,zind]
    n1 <- length(xind)
    n2 <- length(zind)
  } else {
    D <- x@D[,xind,yind,slice]
    mask <- x@mask[xind,yind,slice]
    n1 <- length(xind)
    n2 <- length(yind)
  }
  if(what=="GA"){
  z <- .Fortran("dti2Dga",
                as.double(D),
                as.integer(n1),
                as.integer(n2),
                as.logical(mask),
                fa=double(n1*n2),
                md=double(n1*n2),
                andir=double(3*n1*n2),
                DUPL=FALSE,
                PACKAGE="dti")[c("fa","md","andir")]
  } else {
  z <- .Fortran("dti2Dfa",
                as.double(D),
                as.integer(n1),
                as.integer(n2),
                as.logical(mask),
                fa=double(n1*n2),
                md=double(n1*n2),
                andir=double(3*n1*n2),
                DUPL=FALSE,
                PACKAGE="dti")[c("fa","md","andir")]
   }
   oldpar <- par(mfrow=c(3,3),mar=mar,mgp=mgp,...)
#  now draw information to graphical device
   on.exit(par(oldpar))
   img<-D[1,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   show.image(make.image(65535*img/max(img)))
   title(paste("Dxx: mean",signif(mean(D[mask]),3),"max",signif(max(D[1,,][mask]),3)))
   img<-D[2,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(img))
   title(paste("Dxy: min",signif(min(D[2,,][mask]),3),"max",signif(max(D[2,,][mask]),3)))
   img<-D[3,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(img))
   title(paste("Dxz: min",signif(min(D[3,,][mask]),3),"max",signif(max(D[3,,][mask]),3)))
   show.image(make.image(matrix(z$fa,n1,n2)))
   if(what=="GA"){
   title(paste("Anisotropy Index  (FA)  range:",signif(min(z$fa[mask]),3),"-",
                signif(max(z$fa[mask]),3)))
   } else {
   title(paste("Geodesic Anisotropy (GA)  range:",signif(min(z$fa[mask]),3),"-",
                signif(max(z$fa[mask]),3)))
   }
   img<-D[4,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(65535*img/max(img)))
   title(paste("Dyy: min",signif(min(D[4,,][mask]),3),"max",signif(max(D[4,,][mask]),3)))
   img<-D[5,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(img))
   title(paste("Dyz: min",signif(min(D[5,,][mask]),3),"max",signif(max(D[5,,][mask]),3)))
   andir.image(matrix(z$fa,n1,n2),array(z$andir,c(3,n1,n2)),quant=quant,minanindex=minanindex)
   title(paste("Anisotropy directions"))
   img <- matrix(z$md,n1,n2)
   show.image(make.image(65535*img/max(img)))
   if(what=="GA"){
   title(paste("Mean log diffusivity   range:",signif(min(z$md[mask]),3),"-",
                signif(max(z$md[mask]),3)))
   } else {
   title(paste("Mean diffusivity   range:",signif(min(z$md[mask]),3),"-",
                signif(max(z$md[mask]),3)))
   }
   img<-D[6,,]
   rg<-quantile(img,qrange)
   img[img>rg[2]]<-rg[2]
   img[img<rg[1]]<-rg[1]
   show.image(make.image(65535*img/max(img)))
   title(paste("Dzz: min",signif(min(D[6,,][mask]),3),"max",signif(max(D[6,,][mask]),3)))
   invisible(NULL)
})

setMethod("plot", "dtiData", function(x, y,slice=1, gradient=NULL, view= "axial", show=TRUE,xind=NULL,yind=NULL,zind=NULL, mar=c(3,3,3,.3),mgp=c(2,1,0), ...) {
if(is.null(x@si)) cat("No dwi data yet")
maxsi <- max(x@si)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
if(is.null(gradient)) gradient <- x@s0ind[1]
if(gradient<1||gradient>x@ngrad) {
   warning("gradient number out of range, show s0 image")
   gradient <- x@s0ind[1]
}
adimpro <- require(adimpro)
if (view == "sagittal") {
   if(slice<1||slice>x@ddim[1]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[1]%/%2
   }
   img <- x@si[slice,yind,zind,gradient]
  } else if (view == "coronal") {
   if(slice<1||slice>x@ddim[2]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[2]%/%2
   }
   img <- x@si[xind,slice,zind,gradient]
  } else {
   if(slice<1||slice>x@ddim[3]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[3]%/%2
   }
   img <- x@si[xind,yind,slice,gradient]
  }
  oldpar <- par(mar=mar,mgp=mgp, ...)
  if(adimpro) {
     img <- make.image(65535*img/maxsi)
     if(show) show.image(img,...)
    } else if(show) {
      image(img,...)
    }
    par(oldpar)
    invisible(img)
}
)
setMethod("plot", "dti", function(x, y, ...) cat("No implementation for class dti\n"))

setMethod("plot", "dtiIndices", 
function(x, y, slice=1, view= "axial", method=1, quant=0, minanindex=NULL, show=TRUE, contrast.enh=1,what="FA",xind=NULL,yind=NULL,zind=NULL, mar=c(3,3,3,.3),mgp=c(2,1,0), ...) {
  if(is.null(x@fa)) cat("No anisotropy index yet")
  if(!(method %in% 1:3)) {
      warning("method out of range, reset to 1")
      method <- 1
  }
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  adimpro <- require(adimpro)
  oldpar <- par(mar=mar,mgp=mgp, ...)
  if (view == "sagittal") {
    anindex <- if(what=="GA") x@ga[slice,yind,zind] else x@fa[slice,yind,zind]
    andirection <- x@andir[,slice,yind,zind]
  } else if (view == "coronal") {
    anindex <- if(what=="GA") x@ga[xind,slice,zind] else x@fa[xind,slice,zind]
    andirection <- x@andir[,xind,slice,zind]
  } else {
    anindex <- if(what=="GA") x@ga[xind,yind,slice] else x@fa[xind,yind,slice]
    andirection <- x@andir[,xind,yind,slice]
  }
    anindex[anindex>1]<-0
    anindex[anindex<0]<-0
  if ((method==1) || (method==2)) {
    if(contrast.enh<1&&fa.contrast.enh>0) anindex <- pmin(anindex/contrast.enh,1)
    if(is.null(minanindex)) minanindex <- quantile(anindex,quant,na.rm=TRUE)
    if (diff(range(anindex,na.rm=TRUE)) == 0) minanindex <- 0
    if(method==1) {
      andirection[1,,] <- abs(andirection[1,,])
      andirection[2,,] <- abs(andirection[2,,])
      andirection[3,,] <- abs(andirection[3,,])
    } else {
      ind<-andirection[1,,]<0
      dim(andirection) <- c(3,prod(dim(ind)))
      andirection[,ind] <- - andirection[,ind]
      andirection[2,] <- (1+andirection[2,])/2
      andirection[3,] <- (1+andirection[3,])/2
      dim(andirection) <- c(3,dim(ind))
    }
    andirection <- aperm(andirection,c(2,3,1))
    andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
    if(adimpro) {
      andirection[is.na(andirection)] <- 0
      andirection <- make.image(andirection,gamma=TRUE)
      if(show) show.image(andirection,...)
    } else if(show) {
      dim(anindex) <- dim(andirection)[2:3]
      image(anindex,...)
    }
    par(oldpar)
    invisible(andirection)
  } else if (method==3) {
    if(adimpro) {
      bary[is.na(bary)] <- 0
      bary <- make.image(aperm(bary,c(2,3,1)))
      if(show) show.image(bary,...)
    } else if(show) {
      image(bary[1,,],...)
    }
    par(oldpar)
    invisible(bary)
  } 
})

#
#
#

dtiData <- function(gradient,imagefile,ddim,xind=NULL,yind=NULL,zind=NULL,level=0,mins0value=0,maxvalue=10000,voxelext=c(1,1,1),orientation=c(1,3,5)) {
  if (any(sort((orientation)%/%2) != 0:2)) stop("invalid orientation \n")
  if (dim(gradient)[2]==3) gradient <- t(gradient)
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]
  s0ind <- (1:ngrad)[apply(abs(gradient),2,max)==0] 
  if (!(file.exists(imagefile))) stop("Image file does not exist")
  cat("Start Data reading",date(),proc.time(), "\n")
  zz <- file(imagefile,"rb")
#  si now contains all images (S_0 and S_I), ngrad includes 
#  number of zero gradients

  if (is.null(xind)) xind <- 1:ddim[1]
  if (is.null(yind)) yind <- 1:ddim[2]
  if (is.null(zind)) zind <- 1:ddim[3]
  si <- numeric()
  for (grad in 1:ngrad) {
    sitemp <- readBin(zz,"integer",prod(ddim),2,FALSE)
    dim(sitemp) <- ddim
    si <- c(si,sitemp[xind,yind,zind])
    cat(".")
  }
  close(zz)
  dim(si) <- c(length(xind),length(yind),length(zind),ngrad)
  dimsi <- dim(si)

#  si <- readBin(zz,"integer",prod(ddim)*ngrad,2,FALSE)
#  close(zz)
  cat("Data successfully read",date(),proc.time(), "\n")

#  if (is.null(xind)) xind <- 1:ddim[1]
#  if (is.null(yind)) yind <- 1:ddim[2]
#  if (is.null(zind)) zind <- 1:ddim[3]
#  dim(si) <- c(ddim,ngrad)
#  si <- si[xind,yind,zind,] 
#  dimsi <- dim(si)

#
#   set correct orientation
#
  xyz <- (orientation)%/%2+1
  swap <- orientation-2*(orientation%/%2)
  if(any(xyz!=1:3)) {
      abc <- 1:3
      abc[xyz] <- abc
      si <- aperm(si,c(abc,4))
      swap[xyz] <- swap
      voxelext[xyz] <- voxelext
      dimsi[xyz] <- dimsi[1:3]
      ddim[xyz] <- ddim[1:3]
      gradient[xyz,] <- gradient
  }
  if(swap[1]==1) {
      si <- si[dimsi[1]:1,,,] 
      gradient[1,] <- -gradient[1,]
      }
  if(swap[2]==1) {
      si <- si[,dimsi[2]:1,,]  
      gradient[2,] <- -gradient[2,]
      }
  if(swap[3]==0) {
      si <- si[,,dimsi[3]:1,]    
      gradient[3,] <- -gradient[3,]
      }
#
#   orientation set to radiological convention
#
  si <- .Fortran("initdata",
                 si=as.integer(si),
                 as.integer(dimsi[1]),
                 as.integer(dimsi[2]),
                 as.integer(dimsi[3]),
                 as.integer(dimsi[4]),
                 as.integer(maxvalue),
                 PACKAGE="dti")$si
#  this replaces the content off all voxel with elements <=0 or >maxvalue by 0
     dim(si) <- dimsi
  level <- max(mins0value,level*mean(si[,,,s0ind][si[,,,s0ind]>0])) # set level to level*mean  of positive s_0 values
  ddim0 <- as.integer(ddim)
  ddim <- as.integer(dim(si)[1:3])

  cat("Create auxiliary statistics",date(),proc.time(), " \n")
  btb <- create.designmatrix.dti(gradient)
  rind <- replind(gradient)
  
  invisible(new("dtiData",
                si     = si,
                btb    = btb,
                ngrad  = ngrad, # = dim(btb)[2]
                s0ind  = s0ind, # indices of S_0 images
                replind = rind,
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = xind,
                yind   = yind,
                zind   = zind,
                level  = level,
                voxelext = voxelext,
                orientation = as.integer(c(0,2,5)),
                source = imagefile)
            )
}

readDWIdata <- function(dirlist, format, nslice, gradient, order = NULL,
                        xind=NULL, yind=NULL, zind=NULL,
                        level=0, mins0value=0, maxvalue=10000,
                        voxelext=NULL, orientation=c(1,3,5)) {
  # basic consistency checks
  if (!(format %in% c("DICOM","NIFTI","ANALYZE","AFNI")))
    stop("Cannot handle other formats then DICOM|NIFTI|ANALYZE|AFNI, found:",format)
  if (any(sort((orientation)%/%2) != 0:2)) stop("invalid orientation \n")
  if (dim(gradient)[2]==3) gradient <- t(gradient)
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]
  s0ind <- (1:ngrad)[apply(abs(gradient),2,max)==0] 
  if (is.null(zind)) zind <- 1:nslice

  # generate file list in specified order
  filelist <- NULL
  for (dd in dirlist) filelist <- c(filelist, paste(dd,list.files(dd),sep=""))
  if (format == "DICOM") {
    if (length(filelist) != ngrad * nslice)
      stop("Number of found files does not match ngrad*nslice",length(filelist))
    if (is.null(order)) {
      order <- 1:(ngrad*nslice)
    } else {
      if (length(order) != ngrad*nslice)
        stop("Length of order vector does not match ngrad*nslice")
    }
    dim(order) <- c(nslice,ngrad)
    order <- order[zind,]
    dim(order) <- NULL
    filelist <- filelist[order]
  } else {
    if (length(filelist) != ngrad)
      stop("Number of found files does not match ngrad",length(filelist),"\nPlease provide each gradient cube in a separate file.")
    if (is.null(order)) {
      order <- 1:ngrad
    } else {
      if (length(order) != ngrad)
        stop("Length of order vector does not match ngrad")
    }
    filelist <- filelist[order]
  }
  # read all DICOM files
  cat("Start reading data",date(),proc.time(), "\n")
  si <- numeric()
  cat("\n")
  ddim <- NULL
  for (ff in filelist) {
    cat(".")
    if (format == "DICOM") {
      data <- read.DICOM(ff)
    } else if (format == "NIFTI") {
      data <- read.NIFTI(ff,setmask=FALSE)
    } else if (format == "ANALYZE") {
      data <- read.ANALYZE(ff,setmask=FALSE)
    } else if (format == "AFNI") {
      data <- read.AFNI(ff,setmask=FALSE)
    } 
    if (is.null(ddim)) ddim <- c(data$dim[1:2],nslice,ngrad)
    if (is.null(voxelext)) {
      if (!is.null(data$delta)) {
        if (!prod(voxelext == data$delta))
          warning("Voxel extension",voxelext,"is not found in data:",data$delta)
        voxelext <- data$delta
      } else {
        warning("Voxel extension neither found nor given!")
      }
    }
    if (is.null(xind)) xind <- 1:data$dim[1]
    if (is.null(yind)) yind <- 1:data$dim[2]
    if (format == "DICOM") {
      si <- c(si,extract.data(data)[xind,yind])
    } else {
      ttt <- extract.data(data)
      si <- c(si,ttt[xind,yind,zind,])
    }
  }
  cat("\n")
  dim(si) <- c(length(xind),length(yind),length(zind),ngrad)
  dimsi <- dim(si)
  cat("Data successfully read",date(),proc.time(), "\n")

  # redefine orientation
  xyz <- (orientation)%/%2+1
  swap <- orientation-2*(orientation%/%2)
  if(any(xyz!=1:3)) {
      abc <- 1:3
      abc[xyz] <- abc
      si <- aperm(si,c(abc,4))
      swap[xyz] <- swap
      voxelext[xyz] <- voxelext
      dimsi[xyz] <- dimsi[1:3]
      ddim[xyz] <- ddim[1:3]
      gradient[xyz,] <- gradient
  }
  if(swap[1]==1) {
      si <- si[dimsi[1]:1,,,] 
      gradient[1,] <- -gradient[1,]
      }
  if(swap[2]==1) {
      si <- si[,dimsi[2]:1,,]  
      gradient[2,] <- -gradient[2,]
      }
  if(swap[3]==0) {
      si <- si[,,dimsi[3]:1,]    
      gradient[3,] <- -gradient[3,]
      }
  # orientation set to radiological convention
  si <- .Fortran("initdata",
                 si=as.integer(si),
                 as.integer(dimsi[1]),
                 as.integer(dimsi[2]),
                 as.integer(dimsi[3]),
                 as.integer(dimsi[4]),
                 as.integer(maxvalue),
                 PACKAGE="dti")$si
  # this replaces the content off all voxel with elements <=0 or >maxvalue by 0
  dim(si) <- dimsi
  level <- max(mins0value,level*mean(si[,,,s0ind][si[,,,s0ind]>0])) # set level to level*mean  of positive s_0 values
  ddim0 <- as.integer(ddim)
  ddim <- as.integer(dim(si)[1:3])
    
  cat("Create auxiliary statistics",date(),proc.time(), " \n")
  btb <- create.designmatrix.dti(gradient)
  rind <- replind(gradient)
  
  invisible(new("dtiData",
                si     = si,
                btb    = btb,
                ngrad  = ngrad, # = dim(btb)[2]
                s0ind  = s0ind, # indices of S_0 images
                replind = rind,
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = xind,
                yind   = yind,
                zind   = zind,
                level  = level,
                voxelext = voxelext,
                orientation = as.integer(c(0,2,5)),
                source = paste(dirlist,collapse="|"))
            )
}
#
#
#


dti <- function(object,  ...) cat("This object has class",class(object),"\n")
setGeneric("dti", function(object,  ...) 
standardGeneric("dti"))


dtiTensor <- function(object,  ...) cat("No DTI tensor calculation defined for this class:",class(object),"\n")


setGeneric("dtiTensor", function(object,  ...) standardGeneric("dtiTensor"))

setMethod("dtiTensor","dtiData",
function(object, method="nonlinear",varmethod="replicates",varmodel="local") {
#  available methods are 
#  "linear" - use linearized model (log-transformed)
#  "nonlinear" - use nonlinear model with parametrization according to Koay et.al. (2006)
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  if(method=="linear"){
     ngrad0 <- ngrad - length(s0ind)
     s0 <- object@si[,,,s0ind]
     si <- object@si[,,,-s0ind]
     if(length(s0ind)>1) s0 <- apply(s0,1:3,mean) 
     mask <- s0 > object@level
     mask <- connect.mask(mask)
     dim(s0) <- dim(si) <- NULL
     ttt <- -log(si/s0)
     ttt[is.na(ttt)] <- 0
     ttt[(ttt == Inf)] <- 0
     ttt[(ttt == -Inf)] <- 0
     dim(ttt) <- c(prod(ddim),ngrad0)
     ttt <- t(ttt)
     cat("Data transformation completed ",date(),proc.time(),"\n")

     btbsvd <- svd(object@btb[,-s0ind])
     solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
     D <- solvebtb%*% ttt
     cat("Diffusion tensors generated ",date(),proc.time(),"\n")

     res <- ttt - t(object@btb[,-s0ind]) %*% D
     rss <- res[1,]^2
     for(i in 2:ngrad0) rss <- rss + res[i,]^2
     dim(rss) <- ddim
     sigma2 <- rss/(ngrad0-6)
     D[c(1,4,6),!mask] <- 1e-6
     D[c(2,3,5),!mask] <- 0
     dim(D) <- c(6,ddim)
     dim(res) <- c(ngrad0,ddim)
     cat("Variance estimates generated ",date(),proc.time(),"\n")
     th0 <- array(s0,object@ddim)
     th0[!mask] <- 0
     gc()
  } else {
#  method == "nonlinear" 
     ngrad0 <- ngrad
     si <- aperm(object@si,c(4,1:3))
     s0 <- si[s0ind,,,]
     if(length(s0ind)>1) s0 <- apply(s0,2:4,mean)
     dim(s0) <- ddim
     mask <- s0 > object@level
     mask <- connect.mask(mask)
     cat("start nonlinear regression",date(),proc.time(),"\n")
     z <- .Fortran("nlrdtirg",
                as.integer(si),
                as.integer(ngrad),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.logical(mask),
                as.double(object@btb),
                th0=as.double(s0),
                D=double(6*prod(ddim)),
                as.integer(200),
                as.double(1e-6),
                res=double(ngrad*prod(ddim)),
                rss=double(prod(ddim)),
                PACKAGE="dti",DUP=FALSE)[c("th0","D","res","rss")]
     cat("successfully completed nonlinear regression ",date(),proc.time(),"\n")
     dim(z$th0) <- ddim
     dim(z$D) <- c(6,ddim)
     dim(z$res) <- c(ngrad,ddim)
     dim(z$rss) <- ddim
     df <- sum(table(object@replind)-1)
     res <- z$res
     D <- z$D
     rss <- z$rss
     th0 <- z$th0
     rm(z)
     gc()
     cat("Start variance estimation ",date(),proc.time(),"\n")
     if(df<1||varmethod!="replicates"){
        sigma2 <- rss/(ngrad-7)
     } else {
#
#  We may want something more sophisticated here in case of
#  replicated designs !!!
#
        df <- sum(table(object@replind)-1)
        hmax <- max(1,(125/df)^(1/3))
        z <- replvar(si,object@replind)
#
#   need to correct for underestimtion of variances due to 
#   truncation of si to integer
#   standard deviation is underestimated by about 0.385 for sd > 2
#
        z <- (sqrt(z)+0.385)^2
        dim(z) <- ddim
#  adaptive bw to achive approx. 200 degrees of freedom
           sigma2 <- awsvar(z,shape=df,hmax=pmax(1,(125/df)^(1/3)),mask=mask)
     }
  }
  if(varmodel=="global") sigma2 <- array(median(sigma2[sigma2>0]),dim(sigma2))
     cat("successfully completed variance estimation ",date(),proc.time(),"\n")
  lags <- c(5,5,3)
  scorr <- .Fortran("mcorr",as.double(res),
                   as.logical(mask),
                   as.integer(ddim[1]),
                   as.integer(ddim[2]),
                   as.integer(ddim[3]),
                   as.integer(ngrad0),
                   double(prod(ddim)),
                   double(prod(ddim)),
                   scorr = double(prod(lags)),
                   as.integer(lags[1]),
                   as.integer(lags[2]),
                   as.integer(lags[3]),
                   PACKAGE="dti",DUP=FALSE)$scorr
  dim(scorr) <- lags
  scorr[is.na(scorr)] <- 0
  cat("estimated spatial correlations",date(),proc.time(),"\n")
  cat("first order  correlation in x-direction",signif(scorr[2,1,1],3),"\n")
  cat("first order  correlation in y-direction",signif(scorr[1,2,1],3),"\n")
  cat("first order  correlation in z-direction",signif(scorr[1,1,2],3),"\n")

  scorr[is.na(scorr)] <- 0
  bw <- optim(c(2,2,2),corrrisk,method="L-BFGS-B",lower=c(.2,.2,.2),
  upper=c(3,3,3),lag=lags,data=scorr)$par
  bw[bw <= .25] <- 0
  cat("estimated corresponding bandwidths",date(),proc.time(),"\n")

  invisible(new("dtiTensor",
                D     = D,
                th0   = th0,
                sigma = sigma2,
                scorr = scorr, 
                bw = bw, 
                mask = mask,
                hmax = 1,
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
                source = object@source,
                method = method)
            )
})

#
#
#

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


#
#
#

dtiIndices <- function(object, ...) cat("No DTI indices calculation defined for this class:",class(object),"\n")

setGeneric("dtiIndices", function(object, ...) standardGeneric("dtiIndices"))

setMethod("dtiIndices","dtiTensor",
function(object, which) {
  ddim <- object@ddim

  z <- .Fortran("dtiind3D",
                as.double(object@D),
                as.integer(object@ddim[1]),
                as.integer(object@ddim[2]),
                as.integer(object@ddim[3]),
                as.logical(object@mask),
                fa=double(prod(object@ddim)),
                ga=double(prod(object@ddim)),
                md=double(prod(object@ddim)),
                andir=double(3*prod(object@ddim)),
                bary=double(3*prod(object@ddim)),
                DUPL=FALSE,
                PACKAGE="dti")[c("fa","ga","md","andir","bary")]

  invisible(new("dtiIndices",
                fa = array(z$fa,object@ddim),
                ga = array(z$ga,object@ddim),
                md = array(z$md,object@ddim),
                andir = array(z$andir,c(3,object@ddim)),
                bary = array(z$bary,c(3,object@ddim)),
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = object@ddim,
                ddim0 = object@ddim0,
                voxelext = object@voxelext,
                orientation = object@orientation,
                xind  = object@xind,
                yind  = object@yind,
                zind  = object@zind,
                method = object@method,
                level = object@level,
                source= object@source)
            )
})

clipdti <- function(obj,  ...) cat("Data clipping not defined for this class:",class(object),"\n")


setGeneric("clipdti", function(obj,  ...) standardGeneric("clipdti"))


setMethod("clipdti","dtiData",
function(obj,xind=NULL,yind=NULL,zind=NULL){
if(is.null(xind)) xind <- 1:obj@ddim[1]
if(is.null(yind)) yind <- 1:obj@ddim[2]
if(is.null(zind)) zind <- 1:obj@ddim[3]
invisible(new("dtiData",
                si     = obj@si[xind,yind,zind,],
                btb    = obj@btb,
                ngrad  = obj@ngrad, # = dim(btb)[2]
                s0ind  = obj@s0ind, # indices of S_0 images
                replind = obj@rind,
                ddim   = c(length(xind),length(yind),length(zind)),
                ddim0  = obj@ddim0,
                xind   = obj@xind[xind],
                yind   = obj@yind[yind],
                zind   = obj@zind[zind],
                level  = obj@level,
                voxelext = obj@voxelext,
                orientation = obj@orientation,
                source = obj@imagefile)
            )
})

setMethod("clipdti","dtiTensor",
function(obj,xind=NULL,yind=NULL,zind=NULL){
if(is.null(xind)) xind <- 1:obj@ddim[1]
if(is.null(yind)) yind <- 1:obj@ddim[2]
if(is.null(zind)) zind <- 1:obj@ddim[3]
invisible(new("dtiTensor",
                D     = obj$D[,xind,yind,zind],
                th0   = obj$th0[xind,yind,zind],
                sigma = obj$sigma2[xind,yind,zind],
                scorr = scorr, 
                bw = bw,
                mask = mask[xind,yind,zind],
                hmax = 1,
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                replind = object@replind,
                ddim  = c(length(xind),length(yind),length(zind)),
                ddim0 = object@ddim0,
                xind  = object@xind[xind],
                yind  = object@yind[yind],
                zind  = object@zind[zind],
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                source = object@source,
                method = method)
            )
})

setMethod("clipdti","dtiIndices",function(obj,xind=NULL,yind=NULL,zind=NULL){
if(is.null(xind)) xind <- 1:obj@ddim[1]
if(is.null(yind)) yind <- 1:obj@ddim[2]
if(is.null(zind)) zind <- 1:obj@ddim[3]
  invisible(new("dtiIndices",
                fa = obj$fa[xind,yind,zind],
                ga = obj$ga[xind,yind,zind],
                md = obj$md[xind,yind,zind],
                andir = obj$andir[,xind,yind,zind],
                bary = obj$bary[,xind,yind,zind],
                btb   = object@btb,
                ngrad = object@ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = c(length(xind),length(yind),length(zind)),
                ddim0 = object@ddim0,
                voxelext = object@voxelext,
                orientation = object@orientation,
                xind  = object@xind[xind],
                yind  = object@yind[yind],
                zind  = object@zind[zind],
                method = object@method,
                level = object@level,
                source= object@source)
            )
})

extract <- function(obj,  ...) cat("Data extraction not defined for this class:",class(object),"\n")

setGeneric("extract", function(obj,  ...) standardGeneric("extract"))

setMethod("extract","dtiData",function(obj,what="dtiData",xind=NULL,yind=NULL,zind=NULL){
what <- tolower(what) 
if(! is.character(what)) stop("Argument what needs to be character\n")
if(is.null(xind)) xind <- 1:obj@ddim[1]
if(is.null(yind)) yind <- 1:obj@ddim[2]
if(is.null(zind)) zind <- 1:obj@ddim[3]
z <- list(NULL)
if("btb"%in%what) z$btb <- obj@btb
if("s0"%in%what) z$S0 <- obj@si[xind,yind,zind,obj@s0ind]
if("si"%in%what) z$Si <- obj@si[xind,yind,zind,-obj@s0ind]
if("data"%in%what) z$data <- obj@si[xind,yind,zind,]
invisible(z)
})

setMethod("extract","dtiTensor",function(obj,what="dtiTensor",xind=NULL,yind=NULL,zind=NULL){
what <- tolower(what) 
if(! is.character(what)) stop("Argument what needs to be character\n")
if(is.null(xind)) xind <- 1:obj@ddim[1]
if(is.null(yind)) yind <- 1:obj@ddim[2]
if(is.null(zind)) zind <- 1:obj@ddim[3]
n1 <- length(xind)
n2 <- length(yind)
n3 <- length(zind)
needev <- ("fa"%in%what)||("ga"%in%what)||("md"%in%what)||("evalues"%in%what)
needall <- needev && ("andir"%in%what)
z <- list(NULL)
if(needall){
   erg <- .Fortran("dti3Dall",
                as.double(obj@D[,xind,yind,zind]),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.logical(obj@mask[xind,yind,zind]),
                fa=double(n1*n2*n3),
                ga=double(n1*n2*n3),
                md=double(n1*n2*n3),
                andir=double(3*n1*n2*n3),
                ev=double(3*n1*n2*n3),
                DUPL=FALSE,
                PACKAGE="dti")[c("fa","ga","md","andir","ev")]
if("fa"%in%what) z$fa <- array(erg$fa,c(n1,n2,n3))
if("ga"%in%what) z$ga <- array(erg$ga,c(n1,n2,n3))
if("md"%in%what) z$md <- array(erg$md,c(n1,n2,n3))
if("evalues"%in%what) z$evalues <- array(erg$ev,c(3,n1,n2,n3))
if("andir"%in%what) z$andir <- array(erg$andir,c(3,n1,n2,n3))
} else {
if(needev){
ev <- array(.Fortran("dti3Dev",
               as.double(obj@D[,xind,yind,zind]),
               as.integer(n1),
               as.integer(n2),
               as.integer(n3),
               as.logical(obj@mask[xind,yind,zind]),
               ev=double(3*n1*n2*n3),
               DUPL=FALSE,
               PACKAGE="dti")$ev,c(3,n1,n2,n3))
if("fa"%in%what) {
               dd <- apply(ev^2,2:4,sum)
               md <- (ev[1,,,]+ev[2,,,]+ev[3,,,])/3
               sev <- sweep(ev,2:4,md)
               z$fa <- sqrt(1.5*apply(sev^2,2:4,sum)/dd)
}
if("ga"%in%what) {
               sev <- log(ev)
               md <- (sev[1,,,]+sev[2,,,]+sev[3,,,])/3
               sev <- sweep(sev,2:4,md)
               ga <- sqrt(apply(sev^2,2:4,sum))
               ga[is.na(ga)] <- 0
               z$ga <- ga 
}
if("md"%in%what) z$md <- (ev[1,,,]+ev[2,,,]+ev[3,,,])/3
if("evalues"%in%what) z$evalues <- ev
}
if("andir"%in%what){
    z$andir <- array(.Fortran("dti3Dand",
               as.double(obj@D[,xind,yind,zind]),
               as.integer(n1),
               as.integer(n2),
               as.integer(n3),
               as.logical(obj@mask[xind,yind,zind]),
               andir=double(3*n1*n2*n3),
               DUPL=FALSE,
               PACKAGE="dti")$andir,c(3,n1,n2,n3))
}
}
if("tensor"%in%what) z$tensor <- obj@D[,xind,yind,zind] 
if("s0"%in%what) z$s0 <- obj@th0[xind,yind,zind]
if("mask"%in%what) z$mask <- obj@th0[xind,yind,zind]
invisible(z)
})

setMethod("extract","dtiIndices",function(obj,what,xind=NULL,yind=NULL,zind=NULL){
what <- tolower(what) 
if(! is.character(what)) stop("Argument what needs to be character\n")
if(is.null(xind)) xind <- 1:obj@ddim[1]
if(is.null(yind)) yind <- 1:obj@ddim[2]
if(is.null(zind)) zind <- 1:obj@ddim[3]
n1 <- length(xind)
n2 <- length(yind)
n3 <- length(zind)
z <- list(NULL)
if("fa"%in%what) z$fa <- obj$fa[xind,yind,zind]
if("ga"%in%what) z$ga <- obj$ga[xind,yind,zind]
if("md"%in%what) z$md <- obj$md[xind,yind,zind]
if("andir"%in%what) z$andir <- obj$andir[,xind,yind,zind]
if("bary"%in%what) z$bary <- obj$bary[,xind,yind,zind]
if("bary"%in%what) z$bary <- obj$bary[,xind,yind,zind]
invisible(z)
})

show3d <- function(obj,  ...) cat("3D Visualization not implemented for this class:",class(object),"\n")

setGeneric("show3d", function(obj,  ...) standardGeneric("show3d"))

setMethod("show3d","dtiIndices",function(obj,index="FA",xind=NULL,yind=NULL,zind=NULL,center=NULL,method=1,level=0,bgcolor="black",add=FALSE){
if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
index <- tolower(index) 
if(!(index%in%c("fa","ga"))) stop("index should be either 'FA' or 'GA'\n")
if(is.null(xind)) xind <- 1:obj@ddim[1]
if(is.null(yind)) yind <- 1:obj@ddim[2]
if(is.null(zind)) zind <- 1:obj@ddim[3]
n1 <- length(xind)
n2 <- length(yind)
n3 <- length(zind)
vext <- obj@voxelext
if(is.null(center)) center <- c(median(xind),median(yind),median(zind))
corner <- 2*center*vext
ind <- switch(index,"fa"=obj@fa[xind,yind,zind], "ga"=obj@ga[xind,yind,zind])
ind[ind<level] <- 0
ind <- ind*min(vext)
tmean <- array(0,c(3,n1,n2,n3))
tmean[1,,,] <- xind*vext[1]
tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
andir <- obj@andir[,xind,yind,zind]
if(method==1) {
      andir <- abs(andir)
      dim(andir) <- c(3,n1*n2*n3)
    } else {
      ind<-andir[1,,]<0
      dim(andir) <- c(3,n1*n2*n3)
      andir[,ind] <- - andir[,ind]
      andir[2,] <- (1+andir[2,])/2
      andir[3,] <- (1+andir[3,])/2
    }
colorvalues <- rgb(andir[1,],andir[2,],andir[3,])
dim(andir) <- c(3,n1,n2,n3)
andir <- sweep(andir,2:4,ind,"*")
lcoord <- array(0,c(3,2,n1,n2,n3))
lcoord[,1,,,] <-  andir/2+tmean[,,,]
lcoord[,2,,,] <-  -andir/2+tmean[,,,]
dim(lcoord) <- c(3,2*n1*n2*n3)
lcoord <- cbind(lcoord,corner-rep(.1,3),corner+rep(.1,3))
colorvalues <- c(rbind(colorvalues,colorvalues),bgcolor,bgcolor)
if(!add) {
rgl.open()
rgl.bg(color=bgcolor)
}
rgl.lines(lcoord[1,],lcoord[2,],lcoord[3,],color=colorvalues)
invisible(NULL)
})

setMethod("show3d","dtiTensor",function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,method=1,level=0,scale=.25,
bgcolor="black",add=FALSE,subdivide=2,smooth=TRUE,maxobjects=729,...){
if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
if(is.null(nx)) nx <- obj@ddim[1]
if(is.null(ny)) ny <- obj@ddim[2]
if(is.null(nz)) nz <- obj@ddim[3]
n <- nx*ny*nz
if(is.null(center)) center <- floor(obj@ddim/2)
if(nx*ny*nz>maxobjects) {
warning(paste("size of data cube",n," exceeds maximum of",maxobjects,"\n
        central part of specified cube selected"))
if(nz > maxobjects^(1/3)) n3 <- 1 else n3 <- nz
n1 <- n2 <- floor(sqrt(maxobjects/n3))
} else {
n1 <- nx
n2 <- ny
n3 <- nz
}
xind <- (center[1]-(n1%/%2)):(center[1]+(n1%/%2))
yind <- (center[2]-(n2%/%2)):(center[2]+(n2%/%2))
zind <- (center[3]-(n3%/%2)):(center[3]+(n3%/%2))
xind <- xind[xind>0&xind<=obj@ddim[1]]
yind <- yind[yind>0&yind<=obj@ddim[2]]
zind <- zind[zind>0&zind<=obj@ddim[3]]
n1 <- length(xind)
n2 <- length(yind)
n3 <- length(zind)
n <- n1*n2*n3
if(n==0) stop("Empty cube specified")
cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
"\n yind=",min(yind),":",max(yind),
"\n zind=",min(zind),":",max(zind),"\n")
vext <- obj@voxelext
center <- center*vext
D <- obj@D[,xind,yind,zind]
D <- D/max(D)
dim(D) <- c(6,n)
indpos <- (1:n)[D[1,]*D[4,]*D[6,]>0]
tens <- D[c(1,2,3,2,4,5,3,5,6),indpos]
tmean <- array(0,c(3,n1,n2,n3))
tmean[1,,,] <- xind*vext[1]
tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
dim(tmean) <- c(3,n)
tmean <- tmean[,indpos]
z <- extract(obj,c("andir","fa"),xind,yind,zind)
andir <- z$andir
dim(andir) <- c(3,n1*n2*n3)
andir <- andir[,indpos]
fa <- z$fa[indpos]
n <- length(indpos)
if(method==1) {
      andir <- abs(andir)
    } else {
      ind<-andir[1,]<0
      andir[,ind] <- - andir[,ind]
      andir[2,] <- (1+andir[2,])/2
      andir[3,] <- (1+andir[3,])/2
    }
colorvalues <- rgb(andir[1,],andir[2,],andir[3,])
dim(tens) <- c(3,3,n)
if(level>0){
indpos <- (1:n)[fa>level]
tens <- tens[,,indpos]
tmean <- tmean[,indpos]
colorvalues <- colorvalues[indpos]
fa <- fa[indpos]
n <- length(indpos)
}
if(!add) {
rgl.open()
rgl.bg(color=bgcolor)
}
sphere <- subdivision3d(cube3d(smooth=smooth), subdivide)
norm <- sqrt( sphere$vb[1,]^2 + sphere$vb[2,]^2 + sphere$vb[3,]^2 )
for (i in 1:3) sphere$vb[i,] <- sphere$vb[i,]/norm
sphere$vb[4,] <- 1
sphere$normals <- sphere$vb
tens <- scale*tens
par3d(skipRedraw=TRUE)
cat("Start creating rgl-object\n Progress (out of",n,"):")
for(i in 1:n) {
       if((i%/%100)*100==i) cat(i," ")
       if(i==n) par3d(skipRedraw=FALSE)
       plot3d( ell(sphere, tens[,,i], center=tmean[,i]),
                color=colorvalues[i], alpha=fa[i], add = TRUE, ...)
}
cat("\n")
invisible(NULL)
})

ell <- function (sphere, cov, center = c(0, 0, 0)){
# Adapted from package rgl: 3D visualization device system (OpenGL)
# Authors: Daniel Adler, Duncan Murdoch
  translate3d(rotate3d( sphere, matrix=chol(cov)), center[1], center[2], center[3])
}

