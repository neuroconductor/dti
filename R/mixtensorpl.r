mfunpl <- function(par,siq,grad){
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
mfunplwghts <- function(par,siq,grad){
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
           lev <- c(exp(par[1]),-log(sw))
           mix <- w[o]/sw 
           or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par)
}
dwiMixtensorpl <- function(object, ...) cat("No dwiMixtensorpl calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensorpl", function(object,  ...) standardGeneric("dwiMixtensorpl"))

setMethod("dwiMixtensorpl","dtiData",function(object, maxcomp=2, p=40, maxneighb=7, method="mixtensor", reltol=1e-8, maxit=5000,ngc=100, optmethod="Nelder-Mead"){
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
#  sdcoef <- object@sdcoef
#  if(all(sdcoef==0)) {
#    cat("No parameters for model of error standard deviation found\n estimating these parameters\n You may prefer to run sdpar before calling dtiTensor")
#    sdcoef <- sdpar(object,interactive=FALSE)@sdcoef
#  }
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
  siind <- getsiind(siq,sneighbors,grad,maxcomp,maxc=.866)
  order <- array(maxcomp,ddim)
  lev <- array(as.integer(0),c(2,ddim))
#  logarithmic eigen values
  mix <- array(0,c(maxcomp,ddim))
  orient <- array(0,c(2,maxcomp,ddim))
  sigma2 <- array(0,ddim)
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  for(i1 in 1:n1) for(i2 in 1:n2) for(i3 in 1:n3){
     if(mask[i1,i2,i3]){
     mc0 <- maxcomp
     ord <- mc0+1
     for(j in 1:mc0) orient[,j,i1,i2,i3] <- paroforient(grad[siind[j+1,i1,i2,i3],])
#
#   these are the gradient vectors corresponding to minima in sherical coordinates
#
     par <- numeric(2*mc0+1)
     lev[,i1,i2,i3] <- par[1] <- log(-log(siq[i1,i2,i3,siind[2,i1,i2,i3]])*c(.8))
#
#  these is an initial estimate for the eigen-value parameter
#
     par[rep(2*(1:mc0),rep(2,mc0))+c(0,1)] <- orient[,1:mc0,i1,i2,i3] 
     rss <- Inf
     krit <- Inf
     for(k in mc0:1){
        if(k<ord) {
#
#  otherwise we would reanalyze a model
#
        lpar <- 2*k+1;
        z <-switch(method,
                  "mixtensor" = optim(par[1:(2*k+1)],mfunpl,siq=siq[i1,i2,i3,],grad=grad,
                   method=optmethod,control=list(maxit=maxit,reltol=reltol))
)
        value <- z$value 
        rss <- min(z$value,rss)
        zz <- mfunplwghts(z$par[1:lpar],siq[i1,i2,i3,],grad)
        ord <- zz$ord
        ttt <- value+2*(3*ord+1)/(ngrad-3*maxcomp-1)*rss
        par <- zz$par
#
#     use directions corresponding to largest weights as initial directions
#
        if(ttt < krit) {
           krit <- ttt
           order[i1,i2,i3] <- ord
           lev[,i1,i2,i3] <- zz$lev
           mix[,i1,i2,i3] <- if(ord==maxcomp) zz$mix else c(zz$mix,rep(0,maxcomp-ord))
           orient[,1:ord,i1,i2,i3] <- zz$orient
           if(method=="Jian2") p[i1,i2,i3] <- z$par[length(z$par)]
       }
     }
   }
   sigma2[i1,i2,i3] <- rss/(ngrad0-3*mc0-1)
#   cat("order",order[i1,i2,i3],"\n")
#   cat("error variance",sigma2[i1,i2,i3]*s0[i1,i2,i3]^2,"\n")
#   cat("ev",c(exp(lev[1,i1,i2,i3]),0,0)+exp(lev[2,i1,i2,i3]),"\n")
#   cat("mix",mix[,i1,i2,i3],"\n")
#   cat("orient",orient[,,i1,i2,i3],"\n")
#   if(method=="Jian2") cat("p",p[i1,i2,i3],"\n") 
    if(igc<ngc){
       igc <- igc+1
    } else {
       igc <- 0
      cat("gc",gc(),"\n")
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

