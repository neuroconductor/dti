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
gmfunpl0 <- function(par,siq,grad,pen=1e2){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl0g",
         as.double(par),#par(lpar)
         as.double(siq),#s(n)
         as.double(t(grad)),#g(3,n)
         as.integer(m),
         as.integer(lpar),
         as.integer(ngrad),
         double(3*m),# d(3,m)
         double(ngrad*m),# z(n,m)
         double(m*m),# v(m,m)
         double(m),# w(m)
         double(ngrad*m),# dkgj(n,m)
         double(ngrad*m),# dkgj2(n,m)
         double(ngrad*m),# ddkdphig(n,m)
         double(ngrad*m),# ddkdetag(n,m)
         double(3*m),# dddphi(3,m)
         double(3*m),# dddeta(3,m)
         double(m*m),# dvdth(m,m)
         double(m*m*m),# dvdphi(m,m,m)
         double(m*m*m),# dvdeta(m,m,m)
         double(ngrad*m*lpar),# dzdpars(n,m,lpar)
         double(m*lpar),# dwdpars(m,lpar)
         double(ngrad*m),# zs(n,m)
         double(ngrad*m),# work1(n,m)
         double(ngrad*m),# work2(n,m)
         double(ngrad),# scopy(n)
         as.double(pen),# pen
         dfdpar=double(lpar),# dfdpar(lpar)
         PACKAGE="dti")$dfdpar
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
           if(ord<m){
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
selisample <- function(ngrad,maxcomp,nguess,dgrad,maxc){
saved.seed <- .Random.seed
set.seed(1)
isample <- matrix(sample(ngrad,maxcomp*nguess,replace=TRUE),maxcomp,nguess)
ind <- rep(TRUE,nguess)
if(maxcomp>1){
ind <- .Fortran("selisamp",
                as.integer(isample),
                as.integer(nguess),
                as.integer(maxcomp),
                as.double(dgrad),
                as.integer(dim(dgrad)[1]),
                ind = logical(nguess),
                as.double(maxc),
                DUPL=FALSE,
                PACKAGE="dti")$ind 
#for(i in 1:nguess) for(j in 1:(maxcomp-1)) for(k in (j+1):maxcomp){
#if(dgrad[isample[j,i],isample[k,i]]>maxc) ind[i] <- FALSE
#}
.Random.seed <- saved.seed
}
isample[,ind]
}

paroforient <- function(dir){
  theta <- acos(dir[3])
  sth <- sin(theta)
  phi <- 0
  if(sth<1e-8) {
    theta <- 0
  } else {
    z <- dir[1]/sth
    if(abs(z)>=1) {
      phi <- if(z<0) 0 else pi
    } else {
      phi <- acos(z)*sign(dir[2])
    }
    if(phi < 0) phi <- phi+2*pi
  }
  c(theta, phi)
}

getsiind2 <- function(si,mask,sigma2,grad,vico,th,maxcomp=3,maxc=.866,nguess=100){
# assumes dim(grad) == c(ngrad,3)
# assumes dim(si) == c(n1,n2,n3,ngrad)
# SO removed
ngrad <- dim(grad)[1]
nvico <- dim(vico)[1]
nsi <- dim(si)[4]
dgrad <- matrix(abs(grad%*%t(vico)),ngrad,nvico)
dgrad <- dgrad/max(dgrad)
dgradi <- matrix(abs(vico%*%t(vico)),nvico,nvico)
dgradi <- dgradi/max(dgradi)
nth <- length(th)
isample <- selisample(nvico,maxcomp,nth*nguess,dgradi,maxc)
#
#  eliminate configurations with close directions 
#
if(maxcomp>1) { 
nguess <- trunc(dim(isample)[2]/nth) 
isample <- array(isample[,1:(nguess*nth)],c(maxcomp,nguess,nth))
} 
# this provides configurations of initial estimates with minimum angle between 
# directions > acos(maxc)
cat("using ",nguess,"guesses for initial estimates\n")
nvoxel <- prod(dim(si)[-4])
siind <- .Fortran("getsiin2",
         as.double(aperm(si,c(4,1:3))),
         as.double(sigma2),
         as.integer(nsi),
         as.integer(nvoxel),
         as.integer(maxcomp),
         as.double(dgrad),
         as.integer(nvico),
         as.double(th),
         as.integer(nth),
         double(ngrad*nvico),
         as.integer(isample),
         as.integer(nguess),
         double(nsi),
         double(nsi*(maxcomp+2)),
         siind=integer((maxcomp+2)*nvoxel),
         krit=double(prod(dim(si)[1:3])),
         as.integer(maxcomp+2),
         as.logical(mask),
         PACKAGE="dti")[c("siind","krit")]
failed <- (siind$krit^2/ngrad) > (sigma2-1e-10)
if(any(failed[mask])){
print((siind$krit[mask])[failed[mask]])
print(((1:prod(dim(si)[1:3]))[mask])[failed[mask]])
print(sum(failed[mask]))
}
list(siind=array(siind$siind,c(maxcomp+2,dim(si)[-4])),
     krit=array(siind$krit,dim(si)[-4]))
}

dwiMixtensor <- function(object, ...) cat("No dwiMixtensor calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor", function(object,  ...) standardGeneric("dwiMixtensor"))

setMethod("dwiMixtensor","dtiData",function(object, maxcomp=3,  p=40, method="mixtensor", reltol=1e-6, maxit=5000,ngc=1000, optmethod="BFGS", nguess=25*maxcomp^2,msc="BIC"){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
#  
#
  set.seed(1)
  pen <- 1e2
  theta <- .5
  maxc <- .866
  args <- sys.call(-1)
  args <- c(object@call,args)
  prta <- proc.time()[1]
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(5*(1+3*maxcomp)>ngrad0){
     maxcomp <- max(1,trunc((ngrad0-5)/15))
     cat("Maximal number of components reduced to", maxcomp,"due to insufficient
          number of gradient directions\n")
  }
  cat("Start search outlier detection at",date(),"\n")
#
#  replace physically meaningless S_i by mena S_0 values
#
  z <- .Fortran("outlier",#misc.f 
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
  cat("End search outlier detection at",date(),"\n")
  si <- array(as.integer(z$si),c(ddim,ngrad))
  index <- if(z$lindex>0) z$index[1:z$lindex] else numeric(0)
  rm(z)
  gc()
  cat("Start generating auxiliary objects",date(),"\n")
#
#  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
#
  z <- .Fortran("sweeps0",# mixtens.f
                as.integer(si[,,,-s0ind,drop=FALSE]),
                as.integer(si[,,,s0ind,drop=FALSE]),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.integer(ns0),
                as.integer(ngrad0),
                as.integer(object@level),
                siq=double(prod(ddim[1:3])*ngrad0),
                s0=double(prod(ddim[1:3])),
                vsi=double(prod(ddim[1:3])),
                mask=logical(prod(ddim[1:3])),
                DUPL=FALSE,
                PACKAGE="dti")[c("siq","s0","vsi","mask")]
  rm(si)
  s0 <- array(z$s0,ddim[1:3])
  siq <- array(z$siq,c(ddim[1:3],ngrad0))
  sigma2 <- array(z$vsi,ddim[1:3])
  mask <- array(z$mask,ddim[1:3])
  rm(z)
  gc()
  npar <- 1+3*(0:maxcomp)
#
#   compute penalty for model selection, default BIC
#
  penIC <- switch(msc,"AIC"=2*npar/ngrad0,"BIC"=log(ngrad0)*npar/ngrad0,
                  "AICC"=(1+npar/ngrad0)/(1-(npar+2)/ngrad0),
                  log(ngrad0)*npar/ngrad0)
  cat("End generating auxiliary objects",date(),"\n")
#
#  avoid situations where si's are larger than s0
#
  grad <- t(object@gradient[,-s0ind])
#
# initial estimates for eigenvalues
#
  cat("Start search for initial estimates of eigenvalues at",date(),"\n")
#
#   initialize EV paramters
#
  lev <- array(.Fortran("getev0", # mixtens.f
               as.double(aperm(siq,c(4,1:3))),
               as.integer(ngrad0),
               as.integer(ddim[1]),
               as.integer(ddim[2]),
               as.integer(ddim[3]),
               lev=double(2*prod(ddim)),
               DUPL=FALSE,
               PACKAGE="dti")$lev,c(2,ddim))
  mlev <- median(lev[2,,,][mask])
#
#   get grid for searching EV-parameters
#
  theta <- theta*mlev
  minth <- .4*theta
  maxth <- 1.2*theta
  nth <- 5
  th <- seq(minth,maxth,length=nth)
  cat("using theta=",th,"\n")
#
#   determine initial estimates for orientations 
#
  cat("Start search for initial directions at",date(),"\n")
  data("polyeders")
  vert <- icosa3$vertices
# remove redundant directions
  vind <- rep(TRUE,dim(vert)[2])
  vind[vert[1,]<0] <- FALSE
  vind[vert[1,]==0 & vert[2,] <0] <- FALSE
  vind[vert[1,]==0 & vert[2,] == 0 &vert[3,]<0] <- FALSE
  vert <- vert[,vind]
#
#  compute initial estimates (EV from grid and orientations from icosa3$vertices)
#
  siind <- getsiind2(siq,mask,sigma2,grad,t(vert),th,maxcomp,maxc=maxc,nguess=nguess)
  krit <- siind$krit # sqrt(sum of squared residuals) for initial estimates
  siind <- siind$siind # components 1: model order 2: 
                       # grid index for EV 2+(1:m) index of orientations
  cat("Model orders for initial estimates")
  print(table(siind[1,,,]))
  cat("Frequencies of selected initial thetas")
  print(table(siind[2,,,]))
 cat("End search for initial values at",date(),"\n")
  order <- array(0,ddim)
#  logarithmic eigen values
  mix <- array(0,c(maxcomp,ddim))
  orient <- array(0,c(2,maxcomp,ddim))
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  prt0 <- proc.time()[1]
  prta <- prt0-prta
#
#   loop over voxel in volume
#
  for(i1 in 1:n1) for(i2 in 1:n2) for(i3 in 1:n3){ # begin loop
     if(mask[i1,i2,i3]){ # begin mask
#   only analyze voxel within mask
     mc0 <- maxcomp
     ord <- mc0+1
     for(j in 1:mc0) { 
          iv <- siind[j+2,i1,i2,i3]
          if(iv==0) iv <- j # this should never happen
          orient[,j,i1,i2,i3] <- paroforient(vert[,iv])
     }
#
#   these are the gradient vectors corresponding to minima in spherical coordinates
#
     if(method=="mixtensor"){
     par <- numeric(2*mc0+1)
#  initialize EV-parameter
     if(siind[2,i1,i2,i3]>0){
     par[1] <- th[siind[2,i1,i2,i3]]
     } else {
     par[1] <- .001
     }
#   initialize orientations
     par[rep(2*(1:mc0),rep(2,mc0))+c(0,1)] <- orient[,1:mc0,i1,i2,i3]
     } 
     sigmai <- sigma2[i1,i2,i3]
     krit <- log(sigmai)+penIC[1]
#
#  use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
#
     for(k in mc0:1){ # begin model order
        if(k<ord) {
#
#  otherwise we would reanalyze a model
#
        if(method=="mixtensor"){
           lpar <- 2*k+1
#
           if(optmethod=="BFGS"){
                 z <- optim(par[1:(2*k+1)],mfunpl0,gmfunpl0,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method="BFGS",control=list(maxit=maxit,reltol=reltol))
           } else {
              z <- optim(par[1:(2*k+1)],mfunpl0,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method=optmethod,control=list(maxit=maxit,reltol=reltol))
           }
        }         
        value <- z$value 
# thats sum of squared residuals + penalties (w<0 or 0>th or or th > 8)
#
#   estimate of sigma from the best fitting model
#
        if(method=="mixtensor"){
            zz <- mfunplwghts0(z$par[1:lpar],siq[i1,i2,i3,],grad)
        } 
        ord <- zz$ord
#  replace sigmai by best variance estimate from currently best model
        if(any(zz$lev<0)||ord<k){
           ttt <- krit
#   parameters not interpretable reduce order
        } else {
           si2new <- value/(ngrad0-3*ord-1)
           ttt <- log(si2new)+penIC[1+ord]
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
           sigma2[i1,i2,i3] <- si2new
       }
     }
   } # end model order
    if(igc<ngc){
       igc <- igc+1
    } else {
       igc <- 1
       ingc <- ingc+1
       prt1 <- proc.time()[1]
       gc()
       cat("Nr. of voxel",ingc*ngc,"time elapsed:",prta+prt1-prt0,"remaining time:",
            (prt1-prt0)/(ingc*ngc)*(sum(mask)-ingc*ngc),"\n")
    }
  }# end mask
  }# end loop
  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
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
                rotation = object@rotation,
                source = object@source,
                outlier = index,
                scale = 1,
                method = method)
            )
   }
)
