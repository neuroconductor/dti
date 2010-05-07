mfunpl0 <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar+1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl0",as.double(par),
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
gmfunpl0 <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar+1)/2
mm1 <- m-1
ngrad <- dim(grad)[1]
.Fortran("mfunpl0g",as.double(par),#par
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
                double(mm1*ngrad),#qv
                double(2*mm1*ngrad),#dqv
                as.integer(mm1),#m-1
                gradient=double(lpar),#dh
                PACKAGE="dti")$gradient
}
mfunplwghts0 <- function(par,siq,grad){
#
#   get weights for Mixtensor-model and extract parameters
#
lpar <- length(par)
m <- (lpar+1)/2
ngrad <- dim(grad)[1]
w<-.Fortran("mfunpl0",as.double(par),
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
           wm1 <- w[-1]
           o <- order(wm1,decreasing=TRUE)
           ord <- sum(wm1>0)
           if(ord<m-1){
              o <- o[1:ord]
           }
           sw <- sum(w)
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
           mix <- (w[o]/sw)[-1] 
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

setMethod("dwiMixtensorpl0","dtiData",function(object, maxcomp=3, ex=.2,  p=40, maxneighb=7, method="mixtensor", reltol=1e-8, maxit=5000,ngc=100, optmethod="Nelder-Mead"){
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
  cat("ngrad",ngrad,"mean(si)",range(apply(object@si[,,,-object@s0ind],1:3,mean)),"\n")
  z <- .Fortran("outlier",
                as.integer(object@si),
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
  s0 <- si[,,,s0ind]
  if(length(s0ind)>1) s0 <- apply(s0,1:3,mean)
  mask <- s0>0
  si <- si[,,,-s0ind]
  cat("ngrad",ngrad0,"mean(si)",range(apply(si,1:3,mean)),"\n")
  siq <- sweep(si,1:3,s0,"/")
  cat("ngrad",ngrad0,"mean(siq)",range(apply(siq,1:3,mean)),"\n")
  siqmed <- apply(siq,1:3,median)
  siqmed[siqmed<.9] <- .9
  siqmed[siqmed>.99] <- .99
  cat("ngrad",ngrad0,"mean(siq)",range(apply(siq,1:3,mean)),"\n")
  siq <- sweep(siq,1:3,siqmed,pmin)
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
  cat("ngrad",ngrad0,"mean(siq)",range(apply(siq,1:3,mean)),"\n")
  sneighbors <- neighbors(grad,maxneighb)
#
#   determine initial estimates for orientations 
#
  siind <- getsiind(siq,sneighbors,grad,maxcomp,maxc=.866)
#
# initial estimates for eigenvalues
#
  cat("ngrad",ngrad0,"mean(siq)",range(apply(siq,1:3,mean)),"\n")
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
           if(optmethod=="BFGS"){
           zz <- mfunplwghts0(par[1:(2*k+1)],siq=siq[i1,i2,i3,],grad=grad)
           count <- 0
           while(any(zz$mix==0)&count<10) {
              par[1] <- -log(siq[i1,i2,i3,siind[2,i1,i2,i3]])*runif(1)
              cat("k1",k,"par",par,"\n")
              zz <- mfunplwghts0(par[1:(2*k+1)],siq=siq[i1,i2,i3,],grad=grad)
              count <- count+1
              if(count>10) cat("problem with initial values in voxel(",c(i1,i2,i3),"\n")
           }
           z <- optim(par[1:(2*k+1)],mfunpl0,gmfunpl0,siq=siq[i1,i2,i3,],grad=grad,
                   method="BFGS",control=list(maxit=maxit,reltol=reltol))
           } else {
           z <- optim(par[1:(2*k+1)],mfunpl0,siq=siq[i1,i2,i3,],grad=grad,
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
           zz <- mfunplwghts0(z$par[1:lpar],siq[i1,i2,i3,],grad)
        } else {
           zz <- mfunpl2wghts(z$par[1:lpar],siq[i1,i2,i3,],grad,pex=p)
        }
        ord <- zz$ord
        if(any(zz$lev<0)){
           ttt <- krit
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
           mix[,i1,i2,i3] <- if(ord==maxcomp) zz$mix else c(zz$mix,rep(0,maxcomp-ord))
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
                as.integer(object@si),
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
  mask <- s0>0
  si <- si[,,,-s0ind,drop=FALSE]
  siq <- sweep(si,1:3,s0,"/")
  siqmed <- apply(siq,1:3,median)
  siqmed[siqmed<.9] <- .9
  siqmed[siqmed>.99] <- .99
  siq <- sweep(siq,1:3,siqmed,pmin)
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
  sneighbors <- neighbors(grad,maxneighb)
#
#   determine initial estimates for orientations 
#
  siind <- getsiind(siq,sneighbors,grad,maxcomp,maxc=.866)
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
           mix[,i1,i2,i3] <- if(ord==maxcomp) zz$mix else c(zz$mix,rep(0,maxcomp-ord))
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
                           as.integer(object@si),
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
             si <- si[,,,-s0ind,drop=FALSE]
             if (length(s0ind)>1) s0 <- apply(s0, 1:3, mean) else dim(s0) <- dim(s0)[1:3]
             # normalized DW data
             siq <- sweep(si,1:3,s0,"/")
             # heuristics to avoid DWI that are larger than s0.
             siqmed <- apply(siq, 1:3, median)
             siqmed[siqmed < .9] <- .9
             siqmed[siqmed > .99] <- .99
             siq <- sweep(siq,1:3,siqmed,pmin)
             # mask for calculation
             mask <- s0 > 0
             grad <- t(object@gradient[,-s0ind])
             sneighbors <- neighbors(grad, maxneighb)
             siind <- getsiind(siq, sneighbors, grad, maxcomp, maxc=.866)
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
                     as.integer(sneighbors),                 # (maxneighb x ngrad) matrix of nearest gradient neighbors
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
                           p           = array(1, c(1,1,1)),
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
