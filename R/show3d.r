################################################################
#                                                              #
# Section for show3d() functions (public)                      #
#                                                              #
################################################################

show3d <- function(obj,  ...) cat("3D Visualization not implemented for this class:",class(obj),"\n")

setGeneric("show3d", function(obj,  ...) standardGeneric("show3d"))

setMethod("show3d","dtiData", function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,scale=.5,bgcolor="black",add=FALSE,maxobjects=729,what="ADC",minalpha=1,power=1,nn=1,normalize=FALSE,box=FALSE,title=FALSE,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(is.null(nx)) nx <- obj@ddim[1]
  if(is.null(ny)) ny <- obj@ddim[2]
  if(is.null(nz)) nz <- obj@ddim[3]
  n <- nx*ny*nz
  if(is.null(center)) center <- floor(obj@ddim/2)
  if(nx*ny*nz>maxobjects) {
  cat("size of data cube",n," exceeds maximum of",maxobjects,"\n")
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
  obj <- obj[xind,yind,zind]
  vext <- obj@voxelext
  tmean <- array(0,c(3,n1,n2,n3))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  dim(tmean) <- c(3,n)
  radii <- extract(obj,"sb")$Si
  s0 <- extract(obj,"S0")$S0
  if(length(dim(s0))==4) s0 <- apply(s0,1:3,mean)
  radii <- sweep(radii,1:3,s0,"/")
  if(what=="ADC") radii <- -log(radii)
  ngrad <- dim(radii)[length(dim(radii))]
  dim(radii) <- c(length(radii)/ngrad,ngrad)
  radii <- t(radii)
  sscale <- scale
  if(what=="colorcoded") sscale <- 1
  if(normalize){
     minradii <- apply(radii,2,min)
     maxradii <- apply(radii,2,max)
     radii <- sweep(radii,2,minradii,"-")
     radii <- sweep(radii,2,maxradii-minradii,"/")^power*sscale
  } else {
     radii <- radii/max(radii)*sscale
  }
  gradient <- obj@gradient[,-obj@s0ind]
  if(!add) {
     rgl.open()
     par3d(...)
     rgl.bg(color=bgcolor)
     }
  if(what=="colorcoded") {
     polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
     ngrad <- dim(gradient)[2]
     n <- dim(radii)[2]
     cat("radii",dim(radii))
     radii <- matrix(.Fortran("datinter",
                        as.double(radii),
                        as.integer(n),
                        as.double(gradient),
                        as.integer(ngrad),
                        as.double(polyeder$vertices),
                        as.integer(polyeder$nv),
                        as.integer(nn),#number of nearest neighbors
                        double(nn),#auxiliary 
                        integer(nn),#auxiliary 
                        polyradii=double(polyeder$nv*n),
                        DUPL=FALSE,
                        PACKAGE="dti")$polyradii,polyeder$nv,n)
     cat("newradii",dim(radii))
     show3d.cdata(radii,polyeder,centers=tmean,minalpha=minalpha,scale=scale,...)
  } else {
     show3d.data(radii,gradient,centers=tmean,minalpha=minalpha,...)
  }
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d(switch(tolower(what),"data"="observed DWI data","adc"="observed ADC"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(tolower(what),"data"="observed diffusion weighted data","adc"="apparent diffusion coefficients from data"),"\n",
  if(normalize) "normalized","\n")
  invisible(rgl.cur())
})

##############

setMethod("show3d","dtiTensor", function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,method=1,level=0,scale=.5,bgcolor="black",add=FALSE,subdivide=2,maxobjects=729,what="ADC",minalpha=.25,normalize=NULL,box=FALSE,title=FALSE,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(!exists("icosa0")) data("polyeders")
  if(subdivide<0||subdivide>4) subdivide <- 3
  if(is.null(nx)) nx <- obj@ddim[1]
  if(is.null(ny)) ny <- obj@ddim[2]
  if(is.null(nz)) nz <- obj@ddim[3]
  n <- nx*ny*nz
  if(is.null(center)) center <- floor(obj@ddim/2)
  if(nx*ny*nz>maxobjects) {
  cat("size of data cube",n," exceeds maximum of",maxobjects,"\n")
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
  obj <- obj[xind,yind,zind]
  vext <- obj@voxelext
  center <- center*vext
  D <- obj@D
  D <- D/max(D)
  dim(D) <- c(6,n)
  indpos <- (1:n)[D[1,]*D[4,]*D[6,]>0]
  tens <- D[,indpos]
  tmean <- array(0,c(3,n1,n2,n3))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  dim(tmean) <- c(3,n)
  tmean <- tmean[,indpos]
  z <- extract(obj,what=c("andir","fa"))
  maxev <- extract(obj,what="evalues")$evalues[3,,,]
  maxev <- maxev[indpos]
  andir <- z$andir
  dim(andir) <- c(3,n1*n2*n3)
  andir <- andir[,indpos]
  fa <- z$fa[indpos]
  mask <- obj@mask[indpos]
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
  dim(tens) <- c(6,n)
  if(level>0){
    indpos <- (1:n)[(fa>level)&mask]
    tens <- tens[,indpos]
    tmean <- tmean[,indpos]
    colorvalues <- colorvalues[indpos]
    fa <- fa[indpos]
    maxev <- maxev[indpos]
    n <- length(indpos)
  }
  if(is.null(normalize)) normalize <- switch(tolower(what),"tensor"=FALSE,"adc"=TRUE)
  polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  radii <- .Fortran(switch(tolower(what),"tensor"="ellradii","adcradii"),
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(tens),
                    as.integer(n),
                    as.double(scale/2),
                    radii=double(n*polyeder$nv),
                    DUPL=FALSE,
                    PACKAGE="dti")$radii
  dim(radii) <- c(polyeder$nv,n)
  if(normalize){
     minradii <- apply(radii,2,min)
     maxradii <- apply(radii,2,max)
     radii <- sweep(radii,2,minradii,"-")
     radii <- sweep(radii,2,maxradii-minradii,"/")*scale
  } else {
  radii <- radii/max(radii)*scale
  }
  if(!add) {
     rgl.open()
     par3d(...)
     rgl.bg(color=bgcolor)
     }
  show3d.tens(radii,polyeder,centers=tmean,colors=colorvalues,alpha=minalpha+(1-minalpha)*fa)
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d(switch(tolower(what),"tensor"="estimated tensors","adc"="estimated ADC (tensor)"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(tolower(what),"tensor"="estimated tensors","adc"="apparent diffusion coefficients from estimated tensors"),"\n",
  if(obj@hmax>1) paste("smoothed with hmax=",obj@hmax),if(normalize) "normalized","\n")
  invisible(rgl.cur())
})

##############

setMethod("show3d","dtiIndices",function(obj, index="FA", nx=NULL, ny=NULL, nz=NULL, center=NULL, method=1, level=0, bgcolor="black", add=FALSE, lwd=1,box=FALSE,title=FALSE,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  index <- tolower(index) 
  if(!(index%in%c("fa","ga"))) stop("index should be either 'FA' or 'GA'\n")
  if(is.null(center)) center <- floor(obj@ddim/2)
  if(is.null(nx)) nx <- obj@ddim[1]
  if(is.null(ny)) ny <- obj@ddim[2]
  if(is.null(nz)) nz <- obj@ddim[3]
  xind <- (center[1]-(nx%/%2)):(center[1]+(nx%/%2))
  yind <- (center[2]-(ny%/%2)):(center[2]+(ny%/%2))
  zind <- (center[3]-(nz%/%2)):(center[3]+(nz%/%2))
  xind <- xind[xind>0&xind<=obj@ddim[1]]
  yind <- yind[yind>0&yind<=obj@ddim[2]]
  zind <- zind[zind>0&zind<=obj@ddim[3]]
  n1 <- length(xind)
  n2 <- length(yind)
  n3 <- length(zind)
  n <- n1*n2*n3
  vext <- obj@voxelext
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
    ind1 <- andir[1,,]<0
    dim(andir) <- c(3,n1*n2*n3)
    andir[,ind1] <- - andir[,ind1]
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
  lcoord <- cbind(lcoord)
  colorvalues <- c(rbind(colorvalues,colorvalues))
  if(!add) {
    open3d()
    par3d(...)
    rgl.bg(color=bgcolor)
  }
  rgl.lines(lcoord[1,],lcoord[2,],lcoord[3,],color=colorvalues,size=lwd)
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d("Main directions",color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),"Main directions of diffusion estimated from the tensor model\n\n")
  if(box) bbox3d()
  invisible(rgl.cur())
})

##############

setMethod("show3d","dwiQball", function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,scale=0.5,bgcolor="black",add=FALSE,subdivide=3,maxobjects=729,minalpha=1,power=1,normalize=TRUE,box=FALSE,title=FALSE,...){
  if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(!exists("icosa0")) data("polyeders")
  if(subdivide<0||subdivide>4) subdivide <- 3
  if(is.null(nx)) nx <- obj@ddim[1]
  if(is.null(ny)) ny <- obj@ddim[2]
  if(is.null(nz)) nz <- obj@ddim[3]
  n <- nx*ny*nz
  if(is.null(center)) center <- floor(obj@ddim/2)
  if(nx*ny*nz>maxobjects) {
  cat("size of data cube",n," exceeds maximum of",maxobjects,"\n")
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
  obj <- obj[xind,yind,zind]
  vext <- obj@voxelext
  center <- center*vext
  sphcoef <- obj@sphcoef
  dim(sphcoef) <- c(dim(sphcoef)[1],prod(dim(sphcoef)[-1]))
  tmean <- array(0,c(3,n1,n2,n3))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  dim(tmean) <- c(3,n)
  polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  sphdesign <- design.spheven(obj@order,polyeder$vertices,obj@lambda)$design
  radii <- t(sphdesign)%*%sphcoef
  if(normalize){
     minradii <- apply(radii,2,min)
     maxradii <- apply(radii,2,max)
     radii <- sweep(radii,2,minradii,"-")
     radii <- sweep(radii,2,maxradii-minradii,"/")^power*scale

  } else {
  radii <- radii/max(radii)*scale
  }
  if(!add) {
     rgl.open()
     par3d(...)
     rgl.bg(color=bgcolor)
  }
  show3d.odf(radii,polyeder,centers=tmean,minalpha=minalpha,...)
  if(box) bbox3d()
  if(is.character(title)) {
     title3d(title,color="white",cex=1.5)
  } else {
     if(title) title3d(switch(tolower(obj@what),"ODF"="ODF","wODF"="Weighted ODF","aODF"="alternative ODF","adc"="ADC (Sph. Harmonics)"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(tolower(obj@what),"qball"="Estimated orientation density function (Qball)","adc"="estimated apparent diffusion coefficients (sperical harmonics"),"\n",
  if(normalize) "normalized","\n")
  invisible(rgl.cur())
})

################################################################
#                                                              #
# Section for show3d() functions (misc)                        #
#                                                              #
################################################################

show3d.tens <- function(radii,polyeder,centers=NULL,colors=NULL,alpha=1,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   if(is.null(colors)){
      colors <- heat.colors(1)
   } 
   if(length(colors)!=n){
      nc <- length(colors)
      nnc <- n%/%nc+1
      colors <- rep(colors,nnc)[1:n]
   }
   if(is.null(alpha)){
      alpha <- 1
   } 
   if(length(alpha)!=n){
      nc <- length(alpha)
      nnc <- n%/%nc+1
      alpha <- rep(alpha,nnc)[1:n]
   }
   nv <- polyeder$nv
   ni <- polyeder$ni*3
   colors <- t(matrix(colors,n,ni))
   alpha <- t(matrix(alpha,n,ni))
   vertices <- array(polyeder$vertices,c(3,nv,n))
   indices <- matrix(polyeder$indices,c(ni,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
   vertices <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors,alpha=alpha,...)
}

#############

show3d.data <- function(radii,vertices,centers=NULL,minalpha=1,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   maxradii <- apply(radii,2,max)
   alpha <- minalpha+(1-minalpha)*sweep(radii,2,maxradii,"/")
   colors <- rgb(abs(vertices[1,]),abs(vertices[2,]),abs(vertices[3,]))
   if(length(alpha)!=n){
      nc <- length(alpha)
      nnc <- n%/%nc+1
      alpha <- rep(alpha,nnc)[1:n]
   }
   nv <- dim(vertices)[2]
   indices <- gettriangles(vertices)$triangles
   ni <- length(indices)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   colors <- matrix(colors,nv,n)
   alpha <- matrix(alpha,nv,n)
   vertices <- array(vertices,c(3,nv,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(vertices)[2]*dim(centers)[2]")
   vertices <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors[indices],alpha=alpha[indices],...)
#   rgl.points(vertices[1,],vertices[2,],vertices[3,],
#                 color=colors,alpha=alpha,size=2,...)
}

#############

show3d.cdata <- function(radii,polyeder,centers=NULL,minalpha=1,scale=.5,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   vertices <- polyeder$vertices
   colors <- rainbow(1000,start=0,end=.7)[pmax(1,pmin(1000,as.integer(1000*radii)))]
   alpha <- minalpha+(1-minalpha)*radii
   nv <- polyeder$nv
   ni <- polyeder$ni*3
   colors <- matrix(colors,nv,n)
   alpha <- matrix(alpha,nv,n)
   vertices <- array(polyeder$vertices,c(3,nv,n))
   indices <- matrix(polyeder$indices,c(ni,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
   vertices <- sweep(scale*vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors[indices],alpha=alpha[indices],...)
}

#############

show3d.odf <- function(radii,polyeder,centers=NULL,minalpha=1,...){
   if(is.null(centers)){
      centers <- matrix(0,3,1)
      n <- 1
   } else {
      dcenters <- dim(centers)
      if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
      n <- dcenters[2]
   }
   vertices <- polyeder$vertices
   colors <- rgb(abs(vertices[1,]),abs(vertices[2,]),abs(vertices[3,]))
   alpha <- rep(minalpha,length(radii))
   nv <- polyeder$nv
   ni <- polyeder$ni*3
   colors <- matrix(colors,nv,n)
   alpha <- matrix(alpha,nv,n)
   vertices <- array(polyeder$vertices,c(3,nv,n))
   indices <- matrix(polyeder$indices,c(ni,n))
   if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
   vertices <- sweep(vertices,2:3,radii,"*")
   vertices <- sweep(vertices,c(1,3),centers,"+")
   dim(vertices) <- c(3,nv*n)
   indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
   rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                 color=colors[indices],alpha=alpha[indices],...)
}
