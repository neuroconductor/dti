
dwiQball <- function(object,  ...) cat("No DWI Q-ball calculation defined for this class:",class(object),"\n")

setGeneric("dwiSqrtODF", function(object,  ...) standardGeneric("dwiSqrtODF"))

setMethod("dwiSqrtODF","dtiData",function(object,what="sqrtODF",order=4,forder=1,lambda=0){
  args <- sys.call(-1)
  args <- c(object@call,args)
  if (!(what %in% c("ODF","wODF","aODF","ADC"))) {
      stop("what should specify either ODF, wODF, aODF, or ADC\n")
          }
  ngrad <- object@ngrad
  ddim <- object@ddim
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  sdcoef <- object@sdcoef
  z <- sioutlier(object@si,(1:ngrad)%in%s0ind)
  si <- array(z$si,c(ddim,ngrad))
  index <- z$index
  rm(z)
  gc()

  # prepare data including mask
  ngrad0 <- ngrad - length(s0ind)
  s0 <- si[,,,s0ind]
  si <- si[,,,-s0ind]
  if (ns0>1) {
    dim(s0) <- c(prod(ddim),ns0)
    s0 <- s0 %*% rep(1/ns0,ns0)
    dim(s0) <- ddim
  }
  mask <- s0 > object@level
  mask <- connect.mask(mask)
  bvalues <- object@bvalues
  while((order+1)*(order+2)/2*(forder+1) >=ngrad0) order <- order-2
  cat("Using",length(lord),"sperical harmonics\n")
  L <- Lnml(lambda,forder,order)
  zeta <- getzeta(si,bvalue...)
  Kqzeta <- Kofgrad(grad,zeta,forder,order,bvalue)
  # now switch for different cases
  if (what=="ODF") {
     cat("Data transformation started ",format(Sys.time()),"\n")
     dim(s0) <- dim(si) <- NULL
     si[is.na(si)] <- 0
     si[(si == Inf)] <- 0
     si[(si == -Inf)] <- 0
     dim(si) <- c(prod(ddim),ngrad0)
     si <- t(si)
     cat("Data transformation completed ",format(Sys.time()),"\n")

     # get SH design matrix ...
     z <- design.spheven(order,object@gradient[,-s0ind],lambda)
     # ... and estimate coefficients of SH series of ODF
     # see Descoteaux et al. (2007)
     # include FRT(SH) -> P_l(0)
     sicoef <- z$matrix%*% si
     sphcoef <- plzero(order)%*%sicoef
     cat("Estimated coefficients for ODF (order=",order,") ",format(Sys.time()),"\n")
  } else if (what=="wODF") {
     cat("Data transformation started ",format(Sys.time()),"\n")
     dim(s0) <- dim(si) <- NULL
     si <- si/s0
     si[is.na(si)] <- 0
# Regularization following Aganj et al. (2010) delta=1e-3
     ind1 <- si<0
     ind2 <- (si<1e-3)&(!ind1)
     ind3 <- si>1 
     ind4 <- (si>1-1e-3)&(!ind3)
     si[ind1] <- 5e-4
     si[ind2] <- 5e-4+5e2*si[ind2]^2
     si[ind3] <- 1 - 5e-4
     si[ind4] <- 1 - 5e-4 - 5e2*(1-si[ind4])^2
#     si[si>=1] <- 1-.Machine$double.neg.eps
     si <- log( -log(si))
     si[is.na(si)] <- 0
     si[(si == Inf)] <- 0
     si[(si == -Inf)] <- 0
     dim(si) <- c(prod(ddim),ngrad0)
     si <- t(si)
     cat("Data transformation completed ",format(Sys.time()),"\n")

     # get SH design matrix ...
     z <- design.spheven(order,object@gradient[,-s0ind],lambda)
     # ... and estimate coefficients of SH series of ODF
     # see Aganj et al. (2009)
     # include FRT(SH) -> P_l(0)
     sicoef <- z$matrix%*% si
     plz <- plzero(order)/2/pi
     sphcoef <- plz%*%L%*%sicoef
     coef0 <- sphcoef[1,]
     sphcoef[1,] <- 1/2/sqrt(pi)
     sphcoef[-1,] <- sphcoef[-1,]/8/pi
     cat("Estimated coefficients for wODF (order=",order,") ",format(Sys.time()),"\n")
  } else if (what=="aODF") {
     cat("Data transformation started ",format(Sys.time()),"\n")
     dim(s0) <- dim(si) <- NULL
     si <- si/s0
     si[is.na(si)] <- 0
     si[si>=1] <- 1-.Machine$double.neg.eps
     si <- 1/(-log(si))
     si[is.na(si)] <- 0
     si[(si == Inf)] <- 0
     si[(si == -Inf)] <- 0
     dim(si) <- c(prod(ddim),ngrad0)
     si <- t(si)
     cat("Data transformation completed ",format(Sys.time()),"\n")

     # get SH design matrix ...
     z <- design.spheven(order,object@gradient[,-s0ind],lambda)
     # ... and estimate coefficients of SH series of ODF
     # see Descoteaux et al. (2007)
     # include FRT(SH) -> P_l(0)
     sicoef <- z$matrix%*% si
     sphcoef <- plzero(order)%*%sicoef
     cat("Estimated coefficients for aODF (order=",order,") ",format(Sys.time()),"\n")
  } else { #  what == "ADC" 
     cat("Data transformation started ",format(Sys.time()),"\n")
     dim(s0) <- dim(si) <- NULL
     si <- si/s0
     si[is.na(si)] <- 0
     si[si>=1] <- 1-.Machine$double.neg.eps
     si <- -log(si)
#     si <- -log(si)
     si[is.na(si)] <- 0
     si[(si == Inf)] <- 0
     si[(si == -Inf)] <- 0
     dim(si) <- c(prod(ddim),ngrad0)
     si <- t(si)
     cat("Data transformation completed ",format(Sys.time()),"\n")

     # get SH design matrix ...
     z <- design.spheven(order,object@gradient[,-s0ind],lambda)
     # ... and estimate coefficients of SH series of ADC
     sphcoef <- sicoef <- z$matrix%*% si
     cat("Estimated coefficients for ADC expansion in spherical harmonics (order=",order,") ",format(Sys.time()),"\n")
  }
  res <- si - t(z$design) %*% sicoef
  rss <- res[1,]^2
  for(i in 2:ngrad0) rss <- rss + res[i,]^2
  sigma2 <- rss/(ngrad0-length(lord))
  if(what %in% c("ODF","aODF")){
     varcoef <- outer(diag(plzero(order))^2*diag(z$matrix%*%t(z$matrix)),sigma2,"*")
  } else if(what=="wODF"){
     varcoef <- outer(diag(plzero(order))^2*diag(L)^2*diag(z$matrix%*%t(z$matrix)),sigma2,"*")
     varcoef[-1,] <- varcoef[-1,]/256/pi^4
     varcoef[1,] <- 0
  } else {
     varcoef <- outer(diag(z$matrix%*%t(z$matrix)),sigma2,"*")
  }
  dim(sigma2) <- ddim
  sphcoef[,!mask] <- 0
  dim(sphcoef) <- dim(varcoef) <- c((order+1)*(order+2)/2,ddim)
  dim(res) <- c(ngrad0,ddim)
  cat("Variance estimates generated ",format(Sys.time()),"\n")
  th0 <- array(s0,object@ddim)
  th0[!mask] <- 0
  gc()
#
#   get spatial correlation
#
  scorr <- function(res,mask,ddim,ngrad0,lags=c(5,5,3),mc.cores=mc.cores)
  invisible(new("dwiQball",
                call  = args,
                order = as.integer(order),
                forder = 0,
                zeta = 1,
                lambda = lambda,
                sphcoef = sphcoef,
                varsphcoef = varcoef,
                th0   = th0,
                sigma = sigma2,
                scorr = scorr$scorr, 
                bw = scorr$bw, 
                mask = mask,
                hmax = 1,
                gradient = object@gradient,
                bvalue = object@bvalue,
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
                rotation = object@rotation,
                source = object@source,
                outlier = index,
                scale = 0.5,
                what = what)
            )
})
Qlm <- function(l1,l2,a,m1,m2,b){
    sqrt((2*l1+1)*(2*l2+1)*(2*a+1)/4/pi)*
         coupling_3j(2*l1,2*l2,2*a,0,0,0)*
         coupling_3j(2*l1,2*l2,2*a,2*m1,2*m2,2*b)
}
kappan <- function(n,zeta){
sqrt(2/zeta^1.5*gamma(1:(n+1))/gamma(0:n+3/2))
}
hnni <- function(n){
lni <- function(n,i) (-1)^i*gamma(n+1.5)/gamma(n-i+1)/gamma(i+1.5)/gamma(i+1)
hi <- array(0,c(n+1,n+1,2*n+1))
for(n1 in 0:n){ 
   ln1i <- lni(n1,0:n1)
   for(n2 in 0:n){
      ln2i <- lni(n2,0:n2)
      for(i in 0:n1) for(j in 0:n2){
         k <- i+j
         hi[n1+1,n2+1,k+1] <- hi[n1+1,n2+1,k+1] + ln1i[i+1]*ln2i[j+1]
         }
      }
   }
hi
}
g1f1 <- function(n,L,zeta,q=1){
g1 <- array(0,c(L+1,2*n+1))
pi2zetaq <- pi^2*zeta*q^2
for(i in 0:(2*n)) for(j in (0:L)) {
   g1[j+1,i+1] <- gamma(j+i+1.5)*hyperg_1F1(i+j+1.5, 2*j+1.5, - pi2zetaq)
}
g1
}
Inna <- function(n,L,zeta,q=1){
Inn <- array(0,c(L+1,n+1,n+1))
hi <- hnni(n)
g1 <- g1f1(n,L,zeta,q=q)
kn <- kappan(n,zeta)
for(a in 0:L) {
   za <- zeta^(0:L+1.5)*pi^(2*L+.5)*q^(2*(0:L))/4/gamma(2*(0:L)+1.5)
   for(i in 0:n) for(j in 0:n){
       for(k in 0:(i+j)) {
          Inn[,i+1,j+1] <- Inn[,i+1,j+1]+ za*hi[i+1,j+1,k+1]*g1[,k+1] 
          }
       }
   }
Inn
}
Yab <- function(order,gradients){
  order <- as.integer(max(0,order))
  if(order%%2==1){
    warning("maximum order needs to be even, increase order by one")
    order <- order+1
  } 

  # calculate spherical angles theta and phi corresponding to the gradients
  n <- dim(gradients)[2]
  theta <- phi <- numeric(n)
  for( i in 1:n){
    angles <- sphcoord(gradients[,i])
    theta[i] <- angles[1]
    phi[i] <-  angles[2]
  }

  # values of SH on specified spherical angles 
  getsphericalharmonicseven(order,theta,phi)# first index spharmon (a,b)-index 
  #second index - gradient
}
Qlmmat <- function(order){
# l in (0:(L%/%2))*2 
order <- as.integer(max(0,order))
if(order%%2==1){
warning("maximum order needs to be even, increase order by one")
order <- order+1
} 
kseq <- seq(0,order,2)
nind <- (order+1)*(order+2)/2
kseq2 <- seq(0,2*order,2)
nind2 <- (2*order+1)*(order+1)
Qmat <- array(0,c(nind,nind,nind2))
i1 <- i2 <- i3 <- 1
for(l1 in kseq){
   m1seq <- seq(-l1,l1,1)
   for(m1 in m1seq){
      i2 <- 1
      for(l2 in kseq){
         m2seq <- seq(-l2,l2,1)      
         for(m2 in m2seq){
            i3 <- 1
            for(a in kseq2){
               bseq <- seq(-a,a,1)
               for(b in bseq){
                  if(m1+m2+b==0) Qmat[i1,i2,i3] <- Qlm(l1,l2,a,m1,m2,b)
#  otherwise the result is zero
                  i3 <- i3+1
                  }
               }               
               i2 <- i2+1
            }
         }
         i1 <- i1+1
      }
   }
Qmat
}

Kofgrad <- function(grad,zeta,N=1,L=2,q=1){
#
#   compute the kernel K(q|\zeta) for all gradient directions
#   according to Jian Cheng "Nonnegative Definite EAP and ODF Estimation ...
# 
require(gsl)
ngrad <- dim(grad)[2]
if(length(q)<ngrad) {
q <- rep(1,ngrad)
}
Np1 <- N+1
Lind <- (L+1)*(L+2)/2
Lind2 <- (2*L+1)*(L+1)
a <- NULL
for(i in 0:(L)) a <- c(a,rep(2*i,4*i+1))
Kofg <- array(0,c(Lind,N+1,Lind,N+1,ngrad))
Yabm <- Yab(2*L,grad)# dim Lind2, ngrad
Qlmm <- Qlmmat(L)# dim Lind,Lind,Lind2
Inn <- Inna(N,L,zeta,q=q)# dim L+1 (even values only) , N+1, N+1
for(i1 in 1:Lind) for(i2 in 1:i1)
for(n1 in 1:Np1) for(n2 in 1:n1)
for(ig in 1:ngrad){ 
z <- 0
for(ab in 1:Lind2) {
z <- z+4*pi*(-1)^(a[ab]/2)*Inn[a[ab]/2+1,n1,n2]*Qlmm[i1,i2,ab]*Yabm[ab,ig]
}
Kofg[i1,n1,i2,n2,ig] <- z
# use symmetry
Kofg[i2,n1,i1,n2,ig] <- z
Kofg[i2,n2,i1,n1,ig] <- z
Kofg[i1,n2,i2,n1,ig] <- z
}
}
Lnlm <- function(lambda,N,L){
if(length(lambda)!=2) lambda <- rep(lambda[1],2)
lord <- lambda[2]*rep(seq(0,L,2),2*seq(0,L,2)+1)
nord <- lambda[1]*((0:N)*(1:(N+1)))^2     
Lnlm <- array(rep(lord,N+1),c((L+1)*(L+2)/2,N+1))
for(i in 1:(N+1)) Lnlm[,i] <- Lnlm[,i]+nord[i]
Lnlm
}