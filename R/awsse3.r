# This file contains the implementation of dti.smooth() for 
# "dtiData" lines 12--
# and 
# "dtiTensor" lines 225--
# objects. The latter can be done with 
# "Euclidian Metric" lines 237--
# "Riemann Metric" lines 414--

dwi.smooth <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dwi.smooth", function(object, ...) standardGeneric("dwi.smooth"))

setMethod("dwi.smooth", "dtiData", function(object,kstar,lambda=20,sigma2=NULL,dist="SE3",minsb=5,kexp=.4,coils=0,verbose=FALSE){
  args <- sys.call(-1)
  args <- c(object@call,args)
  sdcoef <- object@sdcoef
  if(length(sdcoef)==4||all(sdcoef[5:8]==0)){
  object <- getsdofsb(object,qA0=.1,qA1=.95,nsb=1,level=NULL)
  }
  kappa <- NULL
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
  hseq <- gethseqse3(kstar,grad,kappa,kexp=kexp,vext,dist=dist,verbose=verbose)
  kappa <- hseq$kappa
  nind <- hseq$n
  hseq <- hseq$h
  th <- sb
  if(is.null(sigma2)){
#
  s2inv <- array((object@sdcoef[5]+object@sdcoef[6]*sb)^{-2},dim(sb))
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
    hakt <- hseq[,k]
    kappa0 <- kappa[,k] 
    param <- lkernse3(hakt,kappa0,grad,vext,nind,dist)
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
                as.integer(coils),
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
lkernse3 <- function(h,kappa,grad,vext,n,dist="SE3"){
      ngrad <- dim(grad)[2]
#      n <- min(1000000,prod(1+2*as.integer(h/c(1,vext)))*ngrad^2)
#      if(dist=="SE3"){ nothing else implemented
      if(length(h)<ngrad) h <- rep(h[1],ngrad)
      if(length(kappa)<ngrad) kappa <- rep(kappa[1],ngrad)
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
#      } 
      dim(z$ind) <- c(5,n)
list(h=h,kappa=kappa,grad=grad,ind=z$ind[,1:z$n],w=z$w[1:z$n],nind=z$n)
}


