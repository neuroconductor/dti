################################################################
#                                                              #
# Section for plot() functions                                 #
#                                                              #
################################################################

setMethod("plot", "dwi", function(x, y, ...) cat("No implementation for class dwi\n"))

##############

setMethod("plot", "dtiData", function(x, y,slice=1, gradient=NULL, view= "axial",
                                      show=TRUE, density=FALSE, xind=NULL, yind=NULL,
                                      zind=NULL, mar=par("mar"),mgp=par("mgp"), ...) {
  if(is.null(x@si)) cat("No dwi data yet")
  maxsi <- max(x@si)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])

# reorient plots depending on .dtiopts
  impars <- get(".dtiopts",envir=.dtiOpts)
  rimpars <- adimpro::rimage.options()
  adimpro::rimage.options(swapx=FALSE,swapy=FALSE)
  if(impars$swapx){
     x@si <- x@si[x@ddim[1]:1,,,]
     xind <- sort((x@ddim[1]+1)-xind)
  }
  if(impars$swapy){
     x@si <- x@si[,x@ddim[2]:1,,]
     yind <- sort((x@ddim[2]+1)-yind)
  }
  if(impars$swapz){
     x@si <- x@si[,,x@ddim[3]:1,]
     zind <- sort((x@ddim[3]+1)-zind)
  }
  if(is.null(gradient)) gradient <- x@s0ind[1]
  if(gradient<1||gradient>x@ngrad) {
    warning("gradient number out of range, show s0 image")
    gradient <- x@s0ind[1]
  }
  if(density) {
    z <- density(x@si[xind,yind,zind,gradient])
    if(show) {
      plot(z,main="Density of S0-values")
      lines(c(x@level,x@level),c(0,max(z$y)),col=2)
    }
    return(invisible(z))
  }
  #adimpro <- require(adimpro)
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
  #if(adimpro) {
  img <- img/maxsi
  if(show) adimpro::rimage(img, ...)
  par(oldpar)
  adimpro::rimage.options(swapx=rimpars$swapx,swapy=rimpars$swapy)
  invisible(adimpro::make.image(65535*img))
})

##############

setMethod("plot", "dtiTensor", function(x, y, slice=1, view="axial", quant=0, minfa=NULL,
                                        contrast.enh=1,what="fa", qrange=c(.01,.99),
                                        xind=NULL,yind=NULL,zind=NULL, mar=par("mar"),
                                        mgp=par("mgp"),...) {
  if(is.null(x@D)) cat("No diffusion tensor yet")
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])

  # reorient plots depending on .dtiopts
    impars <- get(".dtiopts",envir=.dtiOpts)
    rimpars <- adimpro::rimage.options()
    adimpro::rimage.options(swapx=FALSE,swapy=FALSE)
  if(impars$swapx){
     x <- x[x@ddim[1]:1,,]
     xind <- sort((x@ddim[1]+1)-xind)
     if(view == "sagittal") slice <- x@ddim[1]+1-slice
  }
  if(impars$swapy){
     x <- x[,x@ddim[2]:1,]
     yind <- sort((x@ddim[2]+1)-yind)
     if(view == "coronal") slice <- x@ddim[2]+1-slice
  }
  if(impars$swapz){
     x <- x[,,x@ddim[3]:1]
     zind <- sort((x@ddim[3]+1)-zind)
     if(view == "axial") slice <- x@ddim[3]+1-slice
  }

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
    z <- .Fortran(C_dti2dga,
                  as.double(D),
                  as.integer(n1*n2),
                  as.integer(mask),
                  fa=double(n1*n2),
                  md=double(n1*n2),
                  andir=double(3*n1*n2))[c("fa","md","andir")]
  } else {
    z <- .Fortran(C_dti2dfa,
                  as.double(D),
                  as.integer(n1*n2),
                  as.integer(mask),
                  fa=double(n1*n2),
                  md=double(n1*n2),
                  andir=double(3*n1*n2))[c("fa","md","andir")]
  }
  oldpar <- par(mfrow=c(3,3),mar=mar,mgp=mgp,...)
  #  now draw information to graphical device
  on.exit(par(oldpar))
  img<-D[1,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  adimpro::rimage(img, ...)
#  adimpro::show.image(adimpro::make.image(65535*img/max(img)))
  title(paste("Dxx: mean",signif(mean(D[mask]),3),"max",signif(max(D[1,,][mask]),3)))
  img<-D[2,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  img[img<rg[1]]<-rg[1]
  img <- (img-rg[1])/(rg[2]-rg[1])
  adimpro::rimage(img, ...)
#  adimpro::show.image(adimpro::make.image(img))
  title(paste("Dxy: min",signif(min(D[2,,][mask]),3),"max",signif(max(D[2,,][mask]),3)))
  img<-D[3,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  img[img<rg[1]]<-rg[1]
  img <- (img-rg[1])/(rg[2]-rg[1])
#  adimpro::show.image(adimpro::make.image(img))
  adimpro::rimage(img, ...)
  title(paste("Dxz: min",signif(min(D[3,,][mask]),3),"max",signif(max(D[3,,][mask]),3)))
  adimpro::rimage(matrix(z$fa,n1,n2), ...)
#  adimpro::show.image(adimpro::make.image(matrix(z$fa,n1,n2)))
  if(what=="GA"){
    title(paste("Anisotropy Index  (GA)  range:",signif(min(z$fa[mask]),3),"-",
                signif(max(z$fa[mask]),3)))
  } else {
    title(paste("Fractional Anisotropy (FA)  range:",signif(min(z$fa[mask]),3),"-",
                signif(max(z$fa[mask]),3)))
  }
  img<-D[4,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  img[img<rg[1]]<-rg[1]
  adimpro::rimage(img, ...)
#  adimpro::show.image(adimpro::make.image(65535*img/max(img)))
  title(paste("Dyy: min",signif(min(D[4,,][mask]),3),"max",signif(max(D[4,,][mask]),3)))
  img<-D[5,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  img[img<rg[1]]<-rg[1]
  img <- (img-rg[1])/(rg[2]-rg[1])
  adimpro::rimage(img, ...)
#  adimpro::show.image(adimpro::make.image(img))
  title(paste("Dyz: min",signif(min(D[5,,][mask]),3),"max",signif(max(D[5,,][mask]),3)))
  andir.image(matrix(z$fa,n1,n2),array(z$andir,c(3,n1,n2)),quant=quant,minfa=minfa)
  title(paste("Anisotropy directions"))
  img <- matrix(z$md,n1,n2)
  adimpro::rimage(img, ...)
#  adimpro::show.image(adimpro::make.image(65535*img/max(img)))
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
  adimpro::rimage(img, ...)
#  adimpro::show.image(adimpro::make.image(65535*img/max(img)))
  title(paste("Dzz: min",signif(min(D[6,,][mask]),3),"max",signif(max(D[6,,][mask]),3)))
  adimpro::rimage.options(swapx=rimpars$swapx,swapy=rimpars$swapy)
  invisible(NULL)
})

##############

setMethod("plot", "dwiMixtensor", function(x, y, slice=1, view="axial", what="fa", minfa=NULL,
                                           identify=FALSE,  xind=NULL,yind=NULL,zind=NULL,
                                           mar=par("mar"),mgp=par("mgp"),...) {
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])

# reorient plots depending on .dtiopts
  impars <- get(".dtiopts",envir=.dtiOpts)
  rimpars <- adimpro::rimage.options()
  adimpro::rimage.options(swapx=FALSE,swapy=FALSE)
  if(impars$swapx){
     x <- x[x@ddim[1]:1,,]
     xind <- sort((x@ddim[1]+1)-xind)
     if(view == "sagittal") slice <- x@ddim[1]+1-slice
  }
  if(impars$swapy){
     x <- x[,x@ddim[2]:1,]
     yind <- sort((x@ddim[2]+1)-yind)
     if(view == "coronal") slice <- x@ddim[2]+1-slice
  }
  if(impars$swapz){
     x <- x[,,x@ddim[3]:1]
     zind <- sort((x@ddim[3]+1)-zind)
     if(view == "axial") slice <- x@ddim[3]+1-slice
  }

  if (view == "sagittal") {
    x <- x[slice,yind,zind]
  } else if (view == "coronal") {
    x <- x[xind,slice,zind]
  } else {
    x <- x[xind,yind,slice]
  }
  what <- tolower(what)
  stats <- extract(x,what)
  oldpar <- par(mfrow=c(1,length(what)),mar=mar,mgp=mgp,...)
  on.exit(par(oldpar))
  if("w0" %in% what){
    img <- drop(stats$w0)
    adimpro::rimage(img, ...)
#    adimpro::show.image(img <- adimpro::make.image(65535*w0))
    title(paste("Isotropic compartment size"))
  }
  if("fa" %in% what){
    img <- drop(stats$fa)
    if(!is.null(minfa)) img[img<minfa] <- 0
    adimpro::rimage(img, ...)
#    adimpro::show.image(img <- adimpro::make.image(65535*fa))
    title(paste("effective FA"))
  }
  if("order" %in% what){
    img <- drop(stats$order)
    img <- img/max(img)
    adimpro::rimage(img, ...)
#    adimpro::show.image(img <- adimpro::make.image(65535*order/max(order)))
    title(paste("Order of mixture (Maximum=",max(img),")"))
  }
  if("eorder" %in% what){
    img <- drop(stats$eorder)
    img <- img/max(img)
    adimpro::rimage(img, ...)
#    adimpro::show.image(img <- adimpro::make.image(65535*eorder/max(eorder)))
    title(paste("Eff. order of mixture (Maximum=",signif(max(img),2),")"))
  }
  if("ev" %in% what){
    img <- drop(stats$ev[1,,,])
    img <- img/max(img)
    adimpro::rimage(img, ...)
#    adimpro::show.image(img <- adimpro::make.image(65535*ev/max(ev)))
    title(paste("Maximal Eigenvalue (Maximum=",signif(max(img),3),")"))
  }
  if(identify){
    xind<-(1:x@ddim[1])
    yind<-(1:x@ddim[2])
    zind<-(1:x@ddim[3])
    adimpro::rimage(1:dim(img)[1],1:dim(img)[2],img, ...)
    adimpro::rimage.options(swapx=rimpars$swapx,swapy=rimpars$swapy)
    identifyFA(view,slice,xind,yind,zind)
  } else {
    par(oldpar)
    adimpro::rimage.options(swapx=rimpars$swapx,swapy=rimpars$swapy)
    invisible(adimpro::make.image(65535*img))
  }
})

##############

setMethod("plot", "dtiIndices", function(x, y, slice=1, view= "axial", method=1, quant=0, minfa=NULL,
                                         show=TRUE, identify=FALSE, density=FALSE, contrast.enh=1,
                                         what="fa", xind=NULL, yind=NULL, zind=NULL,
                                         mar=par("mar"), mgp=par("mgp"), ...) {
  what <- tolower(what)
  if(is.null(x@fa)) cat("No anisotropy index yet")
  if(!(method %in% 1:6)) {
    warning("method out of range, reset to 1")
    method <- 1
  }
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])

# reorient plots depending on .dtiopts
    impars <- get(".dtiopts",envir=.dtiOpts)
  if(impars$swapx){
     x <- x[x@ddim[1]:1,,]
     xind <- sort((x@ddim[1]+1)-xind)
     if(view == "sagittal") slice <- x@ddim[1]+1-slice
  }
  if(impars$swapy){
     x <- x[,x@ddim[2]:1,]
     yind <- sort((x@ddim[2]+1)-yind)
     if(view == "coronal") slice <- x@ddim[2]+1-slice
  }
  if(impars$swapz){
     x <- x[,,x@ddim[3]:1]
     zind <- sort((x@ddim[3]+1)-zind)
     if(view == "axial") slice <- x@ddim[3]+1-slice
  }

  if(density) {
    x <- x[xind,yind,zind]
    z <- density(if(what=="fa") x@fa[x@fa>0] else x@ga[x@ga>0])
    if(show) {
      plot(z,main=paste("Density of positive",what,"-values"))
    }
    return(invisible(z))
  }
  oldpar <- par(mar=mar,mgp=mgp, ...)
  if (view == "sagittal") {
    anindex <- if(what=="ga") tanh(x@ga[slice,yind,zind]) else x@fa[slice,yind,zind]
    if (method == 3) {
      andirection <- x@bary[,slice,yind,zind]
    } else {
      andirection <- x@andir[,slice,yind,zind]
    }
  } else if (view == "coronal") {
    anindex <- if(what=="ga") tanh(x@ga[xind,slice,zind]) else x@fa[xind,slice,zind]
    if (method == 3) {
      andirection <- x@bary[,xind,slice,zind]
    } else {
      andirection <- x@andir[,xind,slice,zind]
    }
  } else {
    anindex <- if(what=="ga") tanh(x@ga[xind,yind,slice]) else x@fa[xind,yind,slice]
    if (method == 3) {
      andirection <- x@bary[,xind,yind,slice]
    } else {
      andirection <- x@andir[,xind,yind,slice]
    }
  }
  anindex[anindex>1] <- 0
  anindex[anindex<0] <- 0
  if ((method==1) || (method==2) || (method==4)) {
    if (contrast.enh < 1 && contrast.enh > 0) anindex <- pmin(anindex/contrast.enh, 1)
    if(is.null(minfa)) minfa <- quantile(anindex,quant,na.rm=TRUE)
    if (diff(range(anindex,na.rm=TRUE)) == 0) minfa <- 0
    if(method==1) {
      andirection[1,,] <- abs(andirection[1,,])
      andirection[2,,] <- abs(andirection[2,,])
      andirection[3,,] <- abs(andirection[3,,])
    } else if (method==2) {
      ind<-andirection[3,,]<0
      dim(andirection) <- c(3,prod(dim(ind)))
      andirection[,ind] <- - andirection[,ind]
      andirection[1,] <- (1+andirection[1,])/2
      andirection[2,] <- (1+andirection[2,])/2
      andirection <- sweep(andirection,2,sqrt(apply(andirection^2,2,sum)),"/")
      dim(andirection) <- c(3,dim(ind))
    } else {
      andirection[1,,] <- andirection[1,,]^2
      andirection[2,,] <- andirection[2,,]^2
      andirection[3,,] <- andirection[3,,]^2
    }
    andirection <- aperm(andirection,c(2,3,1))
    andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minfa)
    #    if(adimpro) {
    andirection[is.na(andirection)] <- 0
    andirection <- adimpro::make.image(andirection,gammatype="ITU")
    if(show) adimpro::show.image(andirection,...)
    if(identify){
      identifyFA(view,slice,xind,yind,zind)
    } else {
      par(oldpar)
      invisible(andirection)
    }
  } else if (method==3) {
    #    if(adimpro) {
    andirection[is.na(andirection)] <- 0
    bary <- adimpro::make.image(aperm(andirection,c(2,3,1)))
    if(show) adimpro::show.image(bary,...)
    if(identify){
      identifyFA(view,slice,xind,yind,zind)
    } else {
      par(oldpar)
      invisible(bary)
    }
  } else if (method==5) {
    #    if(adimpro) {
    andirection[is.na(andirection)] <- 0
    img.hsi.data <- array(0,dim=c(dim(andirection)[2:3],3))
    img.hsi.data[,,1] <- atan2(andirection[2,,],andirection[1,,])
    img.hsi.data[,,1] <- img.hsi.data[,,1] + pi*(img.hsi.data[,,1]<0)
    img.hsi.data[,,2] <- abs(acos(andirection[3,,]))
    img.hsi.data[,,3] <- anindex
    img.hsi <- adimpro::make.image(img.hsi.data,gammatype="ITU",xmode="HSI")
    if(show) adimpro::show.image(img.hsi,...)
    if(identify){
      identifyFA(view,slice,xind,yind,zind)
    } else {
      par(oldpar)
      invisible(img.hsi)
    }
  } else if (method==6) {
    data("colqFA", envir = environment())
    #    if (adimpro) {
    colqFA <- col2rgb(colqFA)/255
    img.data <- array(0, dim=c(dim(anindex), 3))
    for (i in 1:dim(anindex)[1]) { # i dont like for loops in R!
      for (j in 1:dim(anindex)[2]) {
        img.data[i,j,] <- colqFA[, 255 * anindex[i, j] + 1]
      }
    }
    img <- adimpro::make.image(img.data, gammatype="ITU")
    if(show) adimpro::show.image(img, ...)
    if(identify){
      identifyFA(view,slice,xind,yind,zind)
    } else {
      par(oldpar)
      invisible(img)
    }
  }
})

setMethod("plot", "dwiFiber", function(x, y, ...) {
  plot(density(diff(c(x@startind,dim(x@fibers)[1]+1))-1), ...)
  title("Density of fiber lenghts")
  invisible(NULL)
})

setMethod("plot", "dkiIndices", function(x,
                                         y,
                                         slice = 1,
                                         #view = "axial",
                                         what = c("md", "fa", "mk", "mk2",
                                         "k1", "k2", "k3","kaxial", "kradial", "fak"),
                                         #method = 1,
                                         #quant = 0,
                                         #minfa = NULL,
                                         #show = TRUE,
                                         #identify = FALSE,
                                         #density = FALSE,
                                         #contrast.enh = 1,
                                         xind = NULL,
                                         yind = NULL,
                                         #zind = NULL,
                                         mar=par("mar"),
                                         mgp=par("mgp"),
                                         ...) {
  if(is.null(xind)) xind <- 1:x@ddim[1]
  if(is.null(yind)) yind <- 1:x@ddim[2]

  # reorient plots depending on .dtiopts
  impars <- get(".dtiopts",envir=.dtiOpts)
  rimpars <- adimpro::rimage.options()
  adimpro::rimage.options(swapx=FALSE,swapy=FALSE)
  if(impars$swapx){
     x <- x[x@ddim[1]:1,,]
     xind <- sort((x@ddim[1]+1)-xind)
#     if(view == "sagittal") slice <- x@ddim[1]+1-slice
  }
  if(impars$swapy){
     x <- x[,x@ddim[2]:1,]
     yind <- sort((x@ddim[2]+1)-yind)
#     if(view == "coronal") slice <- x@ddim[2]+1-slice
  }
#  if(impars$swapz){
#     x <- x[,,x@ddim[3]:1]
#     zind <- sort((x@ddim[3]+1)-zind)
#     if(view == "axial") slice <- x@ddim[3]+1-slice
#  }

  switch(what,
         md = adimpro::rimage(xind, yind, x@md[xind, yind, slice], ...),
         fa = adimpro::rimage(xind, yind, x@fa[xind, yind, slice], ...),
         mk = adimpro::rimage(xind, yind, x@mk[xind, yind, slice], ...),
         mk2 = adimpro::rimage(xind, yind, x@mk2[xind, yind, slice], ...),
         k1 = adimpro::rimage(xind, yind, x@k1[xind, yind, slice], ...),
         k2 = adimpro::rimage(xind, yind, x@k2[xind, yind, slice], ...),
         k3 = adimpro::rimage(xind, yind, x@k3[xind, yind, slice], ...),
         kaxial = adimpro::rimage(xind, yind, x@kaxial[xind, yind, slice], ...),
         kradial = adimpro::rimage(xind, yind, x@kradial[xind, yind, slice], ...),
         fak = adimpro::rimage(xind, yind, x@fak[xind, yind, slice], ...))
adimpro::rimage.options(swapx=rimpars$swapx,swapy=rimpars$swapy)
invisible(NULL)
})
