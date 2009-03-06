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
setMethod("show3d","dtiTensor", function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,method=1,level=0,scale=1,bgcolor="black",add=FALSE,subdivide=2,maxobjects=729,what="ADC",...){
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
  cat("Visualization of ",n," tensor objects\n")
  polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  radii <- .Fortran(switch(what,"tensor"="ellradii","adcradii"),
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(tens),
                    as.integer(n),
                    as.double(scale/2),
                    radii=double(n*polyeder$nv),
                    DUPL=FALSE,
                    PACKAGE="dti")$radii
  if(!add) {
     rgl.open()
     par3d(...)
     rgl.bg(color=bgcolor)
     }
  show3d.tens(radii,polyeder,centers=tmean,colors=colorvalues,alpha=fa)
  invisible(rgl.cur())
})
setMethod("show3d","dtiData", function(obj,nx=NULL,ny=NULL,nz=NULL,center=NULL,scale=1,bgcolor="black",add=FALSE,maxobjects=729,minalpha=1,...){
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
  ngrad <- dim(radii)[length(dim(radii))]
  dim(radii) <- c(length(radii)/ngrad,ngrad)
  radii <- t(radii)/max(radii)*scale/2
  gradient <- obj@gradient[,-obj@s0ind]
  if(!add) {
     rgl.open()
     par3d(...)
     rgl.bg(color=bgcolor)
     }
  show3d.data(radii,gradient,centers=tmean,colors=colorvalues,minalpha=minalpha,...)
  invisible(rgl.cur())
})
show3d.data <- function(radii,vertices,centers=NULL,thresh=thresh,minalpha=1,...){
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
gettriangles <- function(gradients){
   dgrad <- dim(gradients)
   if(dgrad[2]==3) gradients <- t(gradients)
   ngrad <- dim(gradients)[2]
   ndist <- ngrad*(ngrad-1)/2
   z <- .Fortran("distvert",
                 as.double(gradients),
                 as.integer(ngrad),
                 ab=integer(2*ndist),
                 distab=double(ndist),
                 as.integer(ndist),
                 DUPL=FALSE,
                 PACKAGE="dti")[c("ab","distab")]
   o <- order(z$distab)
   distab <- z$distab[o]
   ab <- matrix(z$ab,2,ndist)[,o]
   z <- .Fortran("triedges",
                 as.integer(ab),
                 as.double(distab),
                 iab=integer(ndist),
                 as.integer(ndist),
                 triangles=integer(3*5*ngrad),
                 ntriangles=as.integer(5*ngrad),
                 DUPL=FALSE,
                 PACKAGE="dti")[c("iab","triangles","ntriangles")]
list(triangles=matrix(z$triangles,3,5*ngrad)[,1:z$ntriangle],
     edges=ab[,z$iab==2])
}
