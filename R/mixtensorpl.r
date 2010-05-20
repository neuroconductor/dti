mfunpl1 <- function(par,siq,grad,pen=1e2){
#
#   evaluate rss for Mixtensor-model (with isotropic component)
#
lpar <- length(par)
m <- (lpar+1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl1",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                as.double(pen),#penalty for negative weights
                double(ngrad*m),#z(ngrad,m) working array
                double(ngrad),#w(ngrad) working array
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
}
gmfunpl1 <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar+1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl1g",as.double(par),#par
                as.double(siq),#siq
                as.double(t(grad)),#grad
                as.integer(m),#m
                as.integer(lpar),#lpar
                as.integer(ngrad),#ngrad
                double(ngrad*m),#z
                double(m),#w
                double(ngrad),#work
                double(1),#erg
                double(ngrad),#fv
                double(lpar*ngrad),#dfv
                double(m*ngrad),#qv
                double(2*m*ngrad),#dqv
                gradient=double(lpar),#dh
                PACKAGE="dti")$gradient
}
mfunplwghts1 <- function(par,siq,grad,pen=1e2){
#
#   get weights for Mixtensor-model (with isotropic component) and extract parameters
#
lpar <- length(par)
m <- (lpar+1)/2
ngrad <- dim(grad)[1]
w<-.Fortran("mfunpl1",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                as.double(pen),#penalty for negative weights
                double(ngrad*m),#z(ngrad,m) working array
                w=double(ngrad),#w(ngrad) working array
                double(1),#residual sum of squares
                PACKAGE="dti")$w[1:m]
           wm1 <- w[-1] # first comonent corresponds to isotropic part
           o <- order(wm1,decreasing=TRUE)
           ord <- sum(wm1>0)
           if(ord<m-1){
              o <- o[1:ord]
           }
           sw <- sum(w[w>0])
           lev <- c(par[1],-log(sw))
           if(ord>0){
           mix <- wm1[o]/sw
           } else {
           mix <- NULL
           } 
           or <- matrix(par[2:lpar],2,m-1)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par)
}
mfunpl0 <- function(par,siq,grad,pen=1e2){
#
#   evaluate rss for Mixtensor-model (without isotropic component)
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl0",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                as.double(pen),#penalty for negative weights
                double(ngrad*m),#z(ngrad,m) working array
                double(ngrad),#w(ngrad) working array
                erg = double(1),#residual sum of squares
                PACKAGE="dti")$erg
}
gmfunpl0 <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl1g",as.double(par),#par
                as.double(siq),#siq
                as.double(t(grad)),#grad
                as.integer(m),#m
                as.integer(lpar),#lpar
                as.integer(ngrad),#ngrad
                double(ngrad*m),#z
                double(m),#w
                double(ngrad),#work
                double(1),#erg
                double(ngrad),#fv
                double(lpar*ngrad),#dfv
                double(m*ngrad),#qv
                double(2*m*ngrad),#dqv
                gradient=double(lpar),#dh
                PACKAGE="dti")$gradient
}
mfunplwghts0 <- function(par,siq,grad,pen=1e2){
#
#   get weights for Mixtensor-model (without isotropic component) and extract parameters 
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
w<-.Fortran("mfunpl0",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                as.double(pen),#penalty for negative weights
                double(ngrad*m),#z(ngrad,m) working array
                w = double(ngrad),#w(ngrad) working array
                double(1),#residual sum of squares
                PACKAGE="dti")$w[1:m]
           o <- order(w,decreasing=TRUE)
           ord <- sum(w>0)
           if(ord<m-1){
              o <- o[1:ord]
           }
           sw <- sum(w[w>0])
           lev <- c(par[1],-log(sw))
           if(ord>0){
           mix <- w[o]/sw
           } else {
           mix <- NULL
           } 
           or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par)
}
mfunpl <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl",as.double(par),
                as.double(siq),
                as.double(t(grad)),
                as.integer(m),
                as.integer(lpar),
                as.integer(ngrad),
                double(ngrad*m),
                double(m),
                double(ngrad),
                erg = double(1),
                PACKAGE="dti")$erg
}
gmfunpl <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunplgr",as.double(par),#par
                as.double(siq),#siq
                as.double(t(grad)),#grad
                as.integer(m),#m
                as.integer(lpar),#lpar
                as.integer(ngrad),#ngrad
                double(ngrad*m),#z
                double(m),#w
                double(ngrad),#work
                double(1),#erg
                double(ngrad),#fv
                double(lpar*ngrad),#dfv
                double(m*ngrad),#qv
                double(2*m*ngrad),#dqv
                gradient=double(lpar),#dh
                PACKAGE="dti")$gradient
}
mfunplwghts <- function(par,siq,grad){
#
#   get weights for Mixtensor-model and extract parameters
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
w<-.Fortran("mfunpl",as.double(par),
                as.double(siq),
                as.double(t(grad)),
                as.integer(m),
                as.integer(lpar),
                as.integer(ngrad),
                double(ngrad*m),
                w=double(m),
                double(ngrad),
                erg = double(1),
                PACKAGE="dti")$w
           o <- order(w,decreasing=TRUE)
           ord <- sum(w>0)
           if(ord<m){
              o <- o[1:ord]
           }
           sw <- sum(w)
           lev <- c(par[1],-log(sw))
           mix <- (w[o]/sw)
           or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par)
}
mfunpl2 <- function(par,siq,grad,pex){
#
#   evaluate rss for Jian-model
#
lpar <- length(par)
m <- (lpar-2)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl2",as.double(par),
                as.double(siq),
                as.double(t(grad)),
                as.integer(m),
                as.double(pex),
                as.integer(lpar),
                as.integer(ngrad),
                double(ngrad*m),
                double(m),
                double(ngrad),
                erg = double(1),
                PACKAGE="dti")$erg
}
mfunpl2wghts <- function(par,siq,grad,pex){
#
#   get weights for Jian-model and extract parameters
#
lpar <- length(par)
m <- (lpar-2)/2
ngrad <- dim(grad)[1]
w<-.Fortran("mfunpl2",as.double(par),
                as.double(siq),
                as.double(t(grad)),
                as.integer(m),
                as.double(pex),
                as.integer(lpar),
                as.integer(ngrad),
                double(ngrad*m),
                w=double(m),
                double(ngrad),
                erg = double(1),
                PACKAGE="dti")$w
           o <- order(w,decreasing=TRUE)
           ord <- sum(w>0)
           if(ord<m){
              o <- o[1:ord]
           }
           lev <- exp(c(par[1],par[2]))
           mix <- if(ord==0) NULL else w[o]
           or <- matrix(par[3:lpar],2,m)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],par[2],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par)
}

dwiMixtensorpl0 <- function(object, ...) cat("No dwiMixtensorpl calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensorpl0", function(object,  ...) standardGeneric("dwiMixtensorpl0"))

setMethod("dwiMixtensorpl0","dtiData",function(object, maxcomp=3, ex=.2,  p=40, maxneighb=7, method="mixtensor", isocomp=TRUE, reltol=1e-8, maxit=5000,ngc=100, optmethod="Nelder-Mead", pen=1e2,maxc=.866){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
  args <- sys.call(-1)
  args <- c(object@call,args)
  prta <- proc.time()[1]
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  z <- .Fortran("outlier",
                as.double(object@si),
                as.integer(prod(ddim)),
                as.integer(ngrad),
                as.logical((1:ngrad)%in%s0ind),
                as.integer(ns0),
                si=integer(prod(ddim)*ngrad),
                index=integer(prod(ddim)),
                lindex=integer(1),
                DUP=FALSE,
                PACKAGE="dti")[c("si","index","lindex")]
  si <- array(z$si,c(ddim,ngrad))
  index <- if(z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  ngrad0 <- ngrad - length(s0ind)
  s0 <- si[,,,s0ind,drop=FALSE]
  if(length(s0ind)>1) s0 <- apply(s0,1:3,mean) else dim(s0) <- dim(s0)[1:3]
  mask <- s0 > object@level
  siq <- si[,,,-s0ind,drop=FALSE]
  dim(siq) <- c(prod(ddim),ngrad-ns0)
  siq[mask,] <- sweep(siq[mask,],1,s0[mask],"/")
  dim(siq) <- c(ddim,ngrad-ns0)
  siqmed <- apply(siq,1:3,median)
  siqmed[siqmed<.9] <- .9
  siqmed[siqmed>.99] <- .99
  siq <- sweep(siq,1:3,siqmed,pmin)
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
#
#   determine initial estimates for orientations 
#
  siind <- getsiind(siq,mask,grad,maxcomp,maxc=maxc)
#
# initial estimates for eigenvalues
#
  lev <- array(.Fortran("getev0",
               as.double(aperm(siq,c(4,1:3))),
               as.integer(ngrad0),
               as.integer(ddim[1]),
               as.integer(ddim[2]),
               as.integer(ddim[3]),
               lev=double(2*prod(ddim)),
               DUPL=FALSE,
               PACKAGE="dti")$lev,c(2,ddim))
  order <- array(0,ddim)
#  logarithmic eigen values
  mix <- array(0,c(maxcomp,ddim))
  orient <- array(0,c(2,maxcomp,ddim))
  sigma2 <- array(0,ddim)
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  prt0 <- proc.time()[1]
  prta <- prt0-prta
  for(i1 in 1:n1) for(i2 in 1:n2) for(i3 in 1:n3){
     if(mask[i1,i2,i3]){
     mc0 <- maxcomp
     ord <- mc0+1
     for(j in 1:mc0) orient[,j,i1,i2,i3] <- paroforient(grad[siind[j+1,i1,i2,i3],])
#
#   these are the gradient vectors corresponding to minima in spherical coordinates
#
     if(method=="mixtensor"){
     par <- numeric(2*mc0+1)
     par[1] <- ex*lev[2,i1,i2,i3]
#
#  these is an initial estimate for the eigen-value parameter
#
     par[rep(2*(1:mc0),rep(2,mc0))+c(0,1)] <- orient[,1:mc0,i1,i2,i3]
     } else {
     par <- numeric(2*mc0+2)
     par[1:2] <- c((1+lev[2,i1,i2,i3]/p),lev[1,i1,i2,i3]/p)
#
#  Initial values for JIAN-Model need to be checked
#
#
#  these is an initial estimate for the eigen-value parameter
#
     par[rep(2*(1:mc0),rep(2,mc0))+c(1,2)] <- orient[,1:mc0,i1,i2,i3]
     }
     sigma2[i1,i2,i3] <- var(siq[i1,i2,i3,])
     krit <- rss <- (ngrad-1)*sigma2[i1,i2,i3]
     krit <- rss <- Inf
     maxcomp0 <- maxcomp
     for(k in mc0:1){
        if(k<ord) {
#
#  otherwise we would reanalyze a model
#
        if(method=="mixtensor"){
           lpar <- 2*k+1
#
#  in case of isocomp==TRUE use estimates from the more restrictive model as
#  initial parameters
#
           if(optmethod=="BFGS"){
              if(!isocomp||k==mc0)
                 z <- optim(par[1:(2*k+1)],mfunpl0,gmfunpl0,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method="BFGS",control=list(maxit=maxit,reltol=reltol))
              if(isocomp) {
                 if(k==mc0) par[1:(2*k+1)] <- z$par[1:(2*k+1)]
                 z <- optim(z$par[1:(2*k+1)],mfunpl1,gmfunpl1,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method="BFGS",control=list(maxit=maxit,reltol=reltol))
              } 
           } else {
#              cat("i1",i1,"i2",i2,"i3",i3,"k",k,"par",par[1:(2*k+1)],"\n")
              if(!isocomp||k==mc0)
              z <- optim(par[1:(2*k+1)],mfunpl0,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method=optmethod,control=list(maxit=maxit,reltol=reltol))
              if(isocomp) {
                 if(k==mc0) par[1:(2*k+1)] <- z$par[1:(2*k+1)]
                 z <- optim(z$par[1:(2*k+1)],mfunpl1,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method=optmethod,control=list(maxit=maxit,reltol=reltol))
              }
           }
        }         
        value <- z$value
        rss <- min(z$value,rss)
        if(method=="mixtensor"){
           if(isocomp) {
              zz <- mfunplwghts1(z$par[1:lpar],siq[i1,i2,i3,],grad)
           } else {
              zz <- mfunplwghts0(z$par[1:lpar],siq[i1,i2,i3,],grad,pen)
           }
        } else {
           zz <- mfunpl2wghts(z$par[1:lpar],siq[i1,i2,i3,],grad,pex=p)
        }
        ord <- zz$ord
#        cat("ord",ord,"krit",krit,"value",value,"rss",rss,"\n")
#        cat("lev",zz$lev,"par",zz$par,"wghts",zz$mix,"\n")
        if(any(zz$lev<0)||ord<k){
           ttt <- krit
           maxcomp0 <- k-1
#   parameters not interpretable reduce order
        } else {
           penalty <- if(isocomp) 2*(3*ord+2)/(ngrad0-3*maxcomp0-2) else 2*(3*ord+1)/(ngrad0-3*maxcomp0-1)
           ttt <- value+penalty*rss
           par <- zz$par
        }
#
#     use directions corresponding to largest weights as initial directions
#
        if(ttt < krit) {
           krit <- ttt
           order[i1,i2,i3] <- ord
           lev[,i1,i2,i3] <- zz$lev
           mix[,i1,i2,i3] <- if(ord==maxcomp) zz$mix else c(zz$mix,rep(0,maxcomp-ord))
           orient[,1:ord,i1,i2,i3] <- zz$orient
           sigma2[i1,i2,i3] <- rss/(ngrad0-3*maxcomp0-1-isocomp)
       }
     }
   }
#   cat("order",order[i1,i2,i3],"\n")
#   cat("error variance",sigma2[i1,i2,i3]*s0[i1,i2,i3]^2,"\n")
#   cat("ev",c(exp(lev[1,i1,i2,i3]),0,0)+exp(lev[2,i1,i2,i3]),"\n")
#   cat("mix",mix[,i1,i2,i3],"\n")
#   cat("orient",orient[,,i1,i2,i3],"\n")
    if(igc<ngc){
       igc <- igc+1
    } else {
       igc <- 0
       ingc <- ingc+1
       prt1 <- proc.time()[1]
       gc()
       cat("Nr. of voxel",ingc*ngc,"time elapsed:",prta+prt1-prt0,"remaining time:",
            (prt1-prt0)/(ingc*ngc)*(sum(mask)-ingc*ngc),"\n")
    }
  }
  }
  invisible(new("dwiMixtensor",
                call   = args,
                ev     = lev,
                mix    = mix,
                orient = orient,
                order  = order,
                p      = p,
                th0    = s0,
                sigma  = sigma2,
                scorr  = array(1,c(1,1,1)), 
                bw     = c(0,0,0), 
                mask   = mask,
                hmax   = 1,
                gradient = object@gradient,
                btb    = object@btb,
                ngrad  = object@ngrad, # = dim(btb)[2]
                s0ind  = object@s0ind,
                replind = object@replind,
                ddim   = object@ddim,
                ddim0  = object@ddim0,
                xind   = object@xind,
                yind   = object@yind,
                zind   = object@zind,
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                source = object@source,
                outlier = index,
                scale = 1,
                method = method)
            )
}
)

dwiMixtensorpl <- function(object, ...) cat("No dwiMixtensorpl calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensorpl", function(object,  ...) standardGeneric("dwiMixtensorpl"))

setMethod("dwiMixtensorpl","dtiData",function(object, maxcomp=3, ex=.2,  p=40, maxneighb=7, method="mixtensor", reltol=1e-8, maxit=5000,ngc=100, optmethod="Nelder-Mead"){
  args <- sys.call(-1)
  args <- c(object@call,args)
  prta <- proc.time()[1]
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  z <- .Fortran("outlier",
                as.double(object@si),
                as.integer(prod(ddim)),
                as.integer(ngrad),
                as.logical((1:ngrad)%in%s0ind),
                as.integer(ns0),
                si=integer(prod(ddim)*ngrad),
                index=integer(prod(ddim)),
                lindex=integer(1),
                DUP=FALSE,
                PACKAGE="dti")[c("si","index","lindex")]
  si <- array(z$si,c(ddim,ngrad))
  index <- if(z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  ngrad0 <- ngrad - length(s0ind)
  s0 <- si[,,,s0ind,drop=FALSE]
  if(length(s0ind)>1) s0 <- apply(s0,1:3,mean) else dim(s0) <- dim(s0)[1:3]
  mask <- s0 > object@level
  siq <- si[,,,-s0ind,drop=FALSE]
  dim(siq) <- c(prod(ddim),ngrad-ns0)
  siq[mask,] <- sweep(siq[mask,],1,s0[mask],"/")
  dim(siq) <- c(ddim,ngrad-ns0)
  siqmed <- apply(siq,1:3,median)
  siqmed[siqmed<.9] <- .9
  siqmed[siqmed>.99] <- .99
  siq <- sweep(siq,1:3,siqmed,pmin)
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
#
#   determine initial estimates for orientations 
#
  siind <- getsiind(siq,mask,grad,maxcomp,maxc=.866)
  order <- array(0,ddim)
  lev <- array(0,c(2,ddim))
#  logarithmic eigen values
  mix <- array(0,c(maxcomp,ddim))
  orient <- array(0,c(2,maxcomp,ddim))
  sigma2 <- array(0,ddim)
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  prt0 <- proc.time()[1]
  prta <- prt0-prta
  siind0 <- siind
  dim(siind0) <- c(maxcomp+1,n1*n2*n3)
  if(any(siind0[2,mask]==0)){
   cat("siindcontains zeros\n")
   return(list(siind=siind,siq=siq,mask=mask,grad=grad,maxcomp=maxcomp))
  }
  for(i1 in 1:n1) for(i2 in 1:n2) for(i3 in 1:n3){
     if(mask[i1,i2,i3]){
     mc0 <- maxcomp
     ord <- mc0+1
     for(j in 1:mc0) {
       if(siind[j+1,i1,i2,i3]<1||siind[j+1,i1,i2,i3]>ngrad0){
        cat("i1",i1,"i2",i2,"i3",i3,"j",j,"mc0",mc0,"\n")
        cat("siind",siind[,i1,i2,i3],"mask",mask[i1,i2,i3],"\n")
       }
       if(is.na(sin(acos(grad[siind[j+1,i1,i2,i3],3])))){
        cat("i1",i1,"i2",i2,"i3",i3,"j",j,"mc0",mc0,"\n")
        cat("siind",siind[j+1,i1,i2,i3],"\n")
       }
       orient[,j,i1,i2,i3] <- paroforient(grad[siind[j+1,i1,i2,i3],])
       }
#
#   these are the gradient vectors corresponding to minima in spherical coordinates
#
     if(method=="mixtensor"){
     par <- numeric(2*mc0+1)
#     lev[,i1,i2,i3] <-  log(-log(siq[i1,i2,i3,siind[2,i1,i2,i3]])*c(.8,.2))
     lev[,i1,i2,i3] <-  -log(siq[i1,i2,i3,siind[2,i1,i2,i3]])*c(ex,1-ex)
     par[1] <- lev[1,i1,i2,i3]
#
#  these is an initial estimate for the eigen-value parameter
#
     par[rep(2*(1:mc0),rep(2,mc0))+c(0,1)] <- orient[,1:mc0,i1,i2,i3]
     } else {
     par <- numeric(2*mc0+2)
     lev[,i1,i2,i3] <- -log(siq[i1,i2,i3,siind[2,i1,i2,i3]])*c(ex,1-ex)
     par[1:2] <- c((1+lev[2,i1,i2,i3]/p),lev[1,i1,i2,i3]/p)
#
#  these is an initial estimate for the eigen-value parameter
#
     par[rep(2*(1:mc0),rep(2,mc0))+c(1,2)] <- orient[,1:mc0,i1,i2,i3] 
     }
     rss <- Inf
     krit <- Inf
     maxcomp0 <- maxcomp
     for(k in mc0:1){
        if(k<ord) {
#
#  otherwise we would reanalyze a model
#
        if(method=="mixtensor"){
           lpar <- 2*k+1
           if(optmethod=="BFGS"){
           zz <- mfunplwghts(par[1:(2*k+1)],siq=siq[i1,i2,i3,],grad=grad)
           count <- 0
#           while(any(zz$lev<0)&count<10) {
           while(any(zz$mix==0)&count<10) {
              par[1] <- -log(siq[i1,i2,i3,siind[2,i1,i2,i3]])*runif(1)
              zz <- mfunplwghts(par[1:(2*k+1)],siq=siq[i1,i2,i3,],grad=grad)
              count <- count+1
              if(count>10) cat("problem with initial values in voxel(",c(i1,i2,i3),"\n")
           }
           z <- optim(par[1:(2*k+1)],mfunpl,gmfunpl,siq=siq[i1,i2,i3,],grad=grad,
                   method="BFGS",control=list(maxit=maxit,reltol=reltol))
           } else {
           z <- optim(par[1:(2*k+1)],mfunpl,siq=siq[i1,i2,i3,],grad=grad,
                   method=optmethod,control=list(maxit=maxit,reltol=reltol))
           }
        } else {
           lpar <- 2*k+2
           z <- optim(par[1:(2*k+2)],mfunpl2,siq=siq[i1,i2,i3,],grad=grad,pex=p,
                   method=optmethod,control=list(maxit=maxit,reltol=reltol))
        }
        value <- z$value^2
        rss <- min(z$value^2,rss)
        if(method=="mixtensor"){
           zz <- mfunplwghts(z$par[1:lpar],siq[i1,i2,i3,],grad)
        } else {
           zz <- mfunpl2wghts(z$par[1:lpar],siq[i1,i2,i3,],grad,pex=p)
        }
        ord <- zz$ord
        if(any(zz$lev<0)){
           ttt <- Inf
           rss <- Inf
           maxcomp0 <- k-1
#   parameters not interpretable reduce order
        } else {
        ttt <- value+2*(3*ord+1)/(ngrad0-3*maxcomp0-1)*rss
        par <- zz$par
        }
#
#     use directions corresponding to largest weights as initial directions
#
        if(ttt < krit) {
           krit <- ttt
           order[i1,i2,i3] <- ord
           lev[,i1,i2,i3] <- zz$lev
           mix[,i1,i2,i3] <- if(length(zz$mix)==maxcomp) zz$mix else c(zz$mix,rep(0,maxcomp-length(zz$mix)))
           orient[,1:ord,i1,i2,i3] <- zz$orient
       }
     }
   }
   sigma2[i1,i2,i3] <- rss/(ngrad0-3*maxcomp0-1)
#   cat("order",order[i1,i2,i3],"\n")
#   cat("error variance",sigma2[i1,i2,i3]*s0[i1,i2,i3]^2,"\n")
#   cat("ev",c(exp(lev[1,i1,i2,i3]),0,0)+exp(lev[2,i1,i2,i3]),"\n")
#   cat("mix",mix[,i1,i2,i3],"\n")
#   cat("orient",orient[,,i1,i2,i3],"\n")
    if(igc<ngc){
       igc <- igc+1
    } else {
       igc <- 0
       ingc <- ingc+1
       prt1 <- proc.time()[1]
       gc()
       cat("Nr. of voxel",ingc*ngc,"time elapsed:",prta+prt1-prt0,"remaining time:",
            (prt1-prt0)/(ingc*ngc)*(sum(mask)-ingc*ngc),"\n")
    }
  }
  }
  invisible(new("dwiMixtensor",
                call   = args,
                ev     = lev,
                mix    = mix,
                orient = orient,
                order  = order,
                p      = p,
                th0    = s0,
                sigma  = sigma2,
                scorr  = array(1,c(1,1,1)), 
                bw     = c(0,0,0), 
                mask   = mask,
                hmax   = 1,
                gradient = object@gradient,
                btb    = object@btb,
                ngrad  = object@ngrad, # = dim(btb)[2]
                s0ind  = object@s0ind,
                replind = object@replind,
                ddim   = object@ddim,
                ddim0  = object@ddim0,
                xind   = object@xind,
                yind   = object@yind,
                zind   = object@zind,
                voxelext = object@voxelext,
                level = object@level,
                orientation = object@orientation,
                source = object@source,
                outlier = index,
                scale = 1,
                method = method)
            )
}
)

dwiMixtensorpl.new <- function(object, ...) cat("No dwiMixtensorpl.new calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensorpl.new", function(object,  ...) standardGeneric("dwiMixtensorpl.new"))

setMethod("dwiMixtensorpl.new",
          "dtiData",
          function(object,
                   maxcomp = 3,
                   maxneighb = 7,
                   reltol = 1e-8,
                   maxit = 5000,
                   optmethod = "Nelder-Mead")
           {
             cat("entering method dwiMixtensor\n")

             if (maxcomp > 4) stop("maximum number of tensor components is 4\n") # otherwise change number of parameters in C code!

             args <- sys.call(-1)
             args <- c(object@call, args)
             ngrad <- object@ngrad
             ddim <- object@ddim
             s0ind <- object@s0ind
             ns0 <- length(s0ind)

             cat("determine outliers ... ")
             z <- .Fortran("outlier",
                           as.double(object@si),
                           as.integer(prod(ddim)),
                           as.integer(ngrad),
                           as.logical((1:ngrad)%in%s0ind),
                           as.integer(ns0),
                           si      = integer(prod(ddim)*ngrad),
                           index   = integer(prod(ddim)),
                           lindex  = integer(1),
                           DUP     = FALSE,
                           PACKAGE = "dti")[c("si","index","lindex")]
             si <- array(z$si, c(ddim, ngrad))
             index <- if (z$lindex>0) z$index[1:z$lindex] else numeric(0)
             rm(z)
             cat("done\n")

             cat("prepare data and initial estimates ... ")
             # prepare data for optim
             s0 <- si[,,,s0ind,drop=FALSE]
             if (length(s0ind)>1) s0 <- apply(s0, 1:3, mean) else dim(s0) <- dim(s0)[1:3]
             # normalized DW data
             mask <- s0 > object@level
             siq <- si[,,,-s0ind,drop=FALSE]
             dim(siq) <- c(prod(ddim),ngrad-ns0)
             siq[mask,] <- sweep(siq[mask,],1,s0[mask],"/")
             dim(siq) <- c(ddim,ngrad-ns0)
             # heuristics to avoid DWI that are larger than s0.
             siqmed <- apply(siq, 1:3, median)
             siqmed[siqmed < .9] <- .9
             siqmed[siqmed > .99] <- .99
             siq <- sweep(siq,1:3,siqmed,pmin)
             # mask for calculation
             mask <- s0 > object@level
             grad <- t(object@gradient[,-s0ind])
             siind <- getsiind(siq, mask, grad, maxcomp, maxc=.866)
             # siind[1,,,] contains number of potential directions 
             # siind[-1,,,] contains indices of grad corresponding to these directions
             cat("done\n")

             cat("optimizing ... ")
             mm <- switch(optmethod, "Nelder-Mead" = 1,
                                     "BFGS"        = 3,
                                     1)
             pl <- 1

             # perform voxelwise optimization and order selection of tensor mixture model
             a <- .C("mixturepl",
                     as.integer(mm),                         # select mixture method
                     as.integer(prod(ddim)),                 # number of voxels
                     as.integer(mask),                       # calculation mask
                     as.double(siq),                         # DWI without s0
                     as.integer(siind),                      # DWI indices of local minima
                     as.integer(ngrad - length(s0ind)),      # number of DWI
                     as.double(grad),                        # gradient directions
                     as.integer(maxcomp),                    # max number of gradient neighbors
                     as.integer(pl),                         # exp for Jian model
                     as.integer(maxit),                      # max number of iterations for optim
                     as.double(reltol),                      # reltol crit for optim
                     order   = integer(prod(ddim)),           # selected order of mixture
                     lev     = double(2*prod(ddim)),         # logarithmic eigenvalues
                     mix     = double(maxcomp*prod(ddim)),   # mixture weights
                     orient  = double(2*maxcomp*prod(ddim)), # phi/theta for all mixture tensors
                     sigma2  = double(prod(ddim)),           # parameter variance ???
                     DUP     = FALSE,
                     PACKAGE = "dti")[c("order", "lev", "mix", "orient", "sigma2")]

             # set dimension attr
             dim(a$order) <- ddim;
             dim(a$lev) <- c(2, ddim);
             dim(a$mix) <- c(maxcomp, ddim);
             dim(a$orient) <- c(2, maxcomp, ddim);
             dim(a$sigma2) <- ddim;
             cat("done\n")

             # create and return new object
             invisible(new("dwiMixtensor",
                           call        = args,
                           ev          = a$lev,
                           mix         = a$mix,
                           orient      = a$orient,
                           order       = a$order,
                           p           = pl,
                           th0         = s0,
                           sigma       = a$sigma2,
                           scorr       = array(1, c(1,1,1)), # ???
                           bw          = c(0,0,0),           # ???
                           mask        = mask,
                           hmax        = 1,                  # ???
                           gradient    = object@gradient,
                           btb         = object@btb,
                           ngrad       = object@ngrad,
                           s0ind       = object@s0ind,
                           replind     = object@replind,
                           ddim        = object@ddim,
                           ddim0       = object@ddim0,
                           xind        = object@xind,
                           yind        = object@yind,
                           zind        = object@zind,
                           voxelext    = object@voxelext,
                           level       = object@level,
                           orientation = object@orientation,
                           source      = object@source,
                           outlier     = index,
                           scale       = 1,
                           method      = "mixtensor")
                       )
           })
