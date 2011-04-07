# This file contains the implementation of dti.smooth() for 
# "dtiData" lines 12--
# and 
# "dtiTensor" lines 225--
# objects. The latter can be done with 
# "Euclidian Metric" lines 237--
# "Riemann Metric" lines 414--

dwi.smooth <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dwi.smooth", function(object, ...) standardGeneric("dwi.smooth"))

setMethod("dwi.smooth", "dtiData", function(object,kstar,lambda=20,sigma2=NULL,dist="SE3",minsb=5,kappa=NULL,coils=0,verbose=FALSE){
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad <- ngrad - ns0
  grad <- object@gradient[,-s0ind]
  sb <- object@si[,,,-s0ind]
  meansb <- apply(sb,1:3,mean)
  mask <- meansb > minsb
  masksb <- array(mask,c(ddim,ngrad))
  vext <- object@voxelext[2:3]/object@voxelext[1]
  hseq <- gethseqse3(kstar,grad,kappa,vext=vext,dist=dist,verbose=verbose)
  th <- sb
  if(is.null(sigma2)){
     s2inv <- array(1,dim(sb))
  } else {
     s2inv <- array(1/sigma2,dim(sb))  
  }
# make it nonrestrictive for the first step
  ni <- s2inv
  if(is.null(kappa))
  kappa <- if(dist=="SE3") 1/sqrt(.25+.17*ngrad) else .28+.17*ngrad
#  kappa <- getkappa(grad)
  prt0 <- Sys.time()
  cat("adaptive smoothing in SE3, kstar=",kstar,if(verbose)"\n" else " ")
  for(k in 1:kstar){
    hakt <- hseq[k]
    kappa0 <- if(dist=="SE3") kappa else kappa/hakt
    param <- lkernse3(hakt,kappa0,grad,vext,dist)
    z <- .Fortran("adasmse3",
                as.double(sb),
                as.double(th),
                ni=as.double(ni),
                as.logical(mask),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.integer(ngrad),
                as.double(lambda),
                as.integer(param$ind),
                as.double(param$w),
                as.integer(param$n),
                thn=double(prod(ddim)*ngrad),
                r=double(prod(ddim)*ngrad),
                as.double(s2inv),
                double(ngrad),
                double(ngrad),
                double(ngrad),
                DUPL=FALSE,
                PACKAGE="dti")[c("ni","thn","r")]
    ni <- z$ni
    th <- z$thn
    r <- z$r
    ind <- masksb&!is.na(r)&ni>s2inv&r<1e4
if(verbose){
   cat("k:",k,"h_k:",hakt,"quartiles of ni",signif(quantile(z$ni/s2inv),3),
  "mean of ni",signif(mean(z$ni/s2inv),3),
  "time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
    cat("quartiles of r",signif(quantile(r[ind]),3),"length of ind",length(ind),"\n")
    } else {
      cat(".")
    }
#    z <- sum(th[ind]/r[ind]*(ni[ind]-s2inv[ind]))/sum(ni[ind]-s2inv[ind])
#    s2inv <- array(1/z/z,dim(sb))
    yy <- th[ind]/r[ind]
     if(is.null(sigma2)){
        xx <- th[ind]
        ww <- ni[ind]-s2inv[ind]
       if(verbose) print(lm(yy~xx,weights=ww))
        rall <- sum(r[ind]*(ni[ind]-s2inv[ind]))/sum(ni[ind]-s2inv[ind])
#    cat("rall",rall,"\n")
        rall <- max(1.92,rall)
        theta <- koayinv(rall,1)
        th2 <- theta^2   
        s2inv[masksb] <- pi/2*((1+th2/2)*
        besselI(th2/4,0,TRUE)+th2/2*besselI(th2/4,1,TRUE))/pmax(1,th)
        s2inv[s2inv>0.1] <- 0.1
#    cat("quartiles of s2inv",signif(quantile(s2inv),3),"\n")
     }
  }
if(!verbose) cat("\n")
  object@si[,,,-s0ind] <- th
  object@call <- args
  object
}
)

koayinv <- function(r,th0,eps=1e-6){
eps <- max(max(1e-8,eps))
if(r<400){
r <- max(sqrt(pi/(4-pi))+eps,r)
thsq <- th0*th0
db <-(2+thsq)*besselI(thsq/4,0,TRUE)+thsq*besselI(thsq/4,1,TRUE)
ksi <- 2+thsq-pi/8*db*db
th1 <- sqrt(ksi*(1+r*r)-2)
while(abs(th0-th1)>eps){
th0 <- th1
thsq <- th0*th0
db <-(2+thsq)*besselI(thsq/4,0,TRUE)+thsq*besselI(thsq/4,1,TRUE)
ksi <- 2+thsq-pi/8*db*db
th1 <- sqrt(ksi*(1+r*r)-2)
#cat(ksi,r,th0,th1,"\n")
}
} else {
th1 <- r*0.9999953
}
th1
}
lkernse3 <- function(h,kappa,grad,vext,dist="SE3"){
      ngrad <- dim(grad)[2]
      n <- min(1000000,prod(1+2*as.integer(h/c(1,vext)))*ngrad^2)
      if(dist=="SE3"){
      z <- .Fortran("lkse3",
                    as.double(h),
                    as.double(kappa),
                    as.double(grad),
                    double(3*ngrad),
                    double(2*ngrad),
                    as.integer(ngrad),
                    as.double(vext),
                    ind=integer(5*n),
                    w=double(n),
                    n=as.integer(n),
                    DUPL=FALSE,
                    PACKAGE="dti")[c("ind","w","n")]
      } else {
      z <- .Fortran("lockse3",
                    as.double(h),
                    as.double(kappa),
                    as.double(grad),
                    as.integer(ngrad),
                    as.double(vext),
                    ind=integer(5*n),
                    w=double(n),
                    n=as.integer(n),
                    DUPL=FALSE,
                    PACKAGE="dti")[c("ind","w","n")]
      }
      dim(z$ind) <- c(5,n)
list(h=h,kappa=kappa,grad=grad,ind=z$ind[,1:z$n],w=z$w[1:z$n],nind=z$n)
}

smoothdwi <- function(object,kstar,lambda=20,dist="SE3",minsb=5){
  args <- sys.call(-1)
  args <- c(object@call,args)
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad <- ngrad - ns0
  grad <- object@gradient[,-s0ind]
  sb <- object@si[,,,-s0ind]
  meansb <- apply(sb,1:3,mean)
  mask <- meansb > minsb
  masksb <- array(mask,c(ddim,ngrad))
  vext <- object@voxelext[2:3]/object@voxelext[1]
  hseq <- gethseqse3(kstar,grad,vext=vext,dist=dist)
  th <- sb
  s2inv <- array(1e-10,dim(sb))
# make it nonrestrictive for the first step
  ni <- s2inv
  kappa <- if(dist=="SE3") 1/sqrt(.25+.17*ngrad) else .28+.17*ngrad
#  kappa <- getkappa(grad)
  prt0 <- Sys.time()
  for(k in 1:kstar){
    hakt <- hseq[k]
    kappa0 <- if(dist=="SE3") kappa else kappa/hakt
    param <- lkernse3(hakt,kappa0,grad,vext,dist)
    z <- .Fortran("adasmse3",
                as.double(sb),
                as.double(th),
                ni=as.double(ni),
                as.logical(mask),
                as.integer(ddim[1]),
                as.integer(ddim[2]),
                as.integer(ddim[3]),
                as.integer(ngrad),
                as.double(lambda),
                as.integer(param$ind),
                as.double(param$w),
                as.integer(param$n),
                thn=double(prod(ddim)*ngrad),
                r=double(prod(ddim)*ngrad),
                as.double(s2inv),
                double(ngrad),
                double(ngrad),
                double(ngrad),
                DUPL=FALSE,
                PACKAGE="dti")[c("ni","thn","r")]
    ni <- z$ni
    th <- z$thn
    r <- z$r
cat("k:",k,"h_k:",hakt,"quartiles of ni",signif(quantile(z$ni/s2inv),3),
  "mean of ni",signif(mean(z$ni/s2inv),3),
  "time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
    ind <- masksb&!is.na(r)&ni>s2inv&r<1e4
    cat("quartiles of r",signif(quantile(r[ind]),3),"length of ind",length(ind),"\n")
    rall <- sum(r[ind]*(ni[ind]-s2inv[ind]))/sum(ni[ind]-s2inv[ind])
    cat("rall",rall,"\n")
    rall <- max(1.92,rall)
    theta <- koayinv(rall,1)
    th2 <- theta^2
    s2inv[masksb] <- pi/2*exp(-th2/4)*
    ((1+th2/2)*besselI(th2/4,0)+th2/2*besselI(th2/4,1))/pmax(1,th)
    s2inv[s2inv>0.1] <- 0.1
    cat("quartiles of s2inv",signif(quantile(s2inv),3),"\n")
  }
  object@si[,,,-s0ind] <- th
  object@call <- args
  object
}

