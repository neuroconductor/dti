pmatrix <- function(x, FUN, ..., mc.cores = setCores(,reprt=FALSE)){
     require(parallel)
     cl <- makeCluster(mc <- mc.cores)
#
#   parCapply does not pass dimension attributes to FUN !!!!
#
     z <- parCapply(cl, x, FUN, ...)
     stopCluster(cl)
     z
}

plmatrix <- function(x, FUN, ..., mc.cores = setCores(,reprt=FALSE)){
dx <- dim(x)[2]
if(mc.cores>dx) mc.cores <- dx
n <- trunc((dx-1)/mc.cores)+1
lx <- list(NULL)
for(i in 1:(mc.cores-1)) lx[[i]] <- x[,(i-1)*n+1:n]
lx[[mc.cores]] <- x[,((mc.cores-1)*n+1):dx]
cl <- makeCluster(mc <- mc.cores)
lz <- parLapply(cl, lx, FUN , ...)
stopCluster(cl)
z <- matrix(0,length(lz[[1]])/n, dx)
for(i in 1:(mc.cores-1)) z[,(i-1)*n+1:n] <- lz[[i]]
z[,((mc.cores-1)*n+1):dx] <- lz[[mc.cores]]
z
}

pnlrdtirg <- function(si,btb,sdcoef,s0ind,ngrad){
   ns0 <- length(s0ind)
   nvox <- length(si)/ngrad
   dim(si) <- c(ngrad,length(si)/ngrad)
   s0 <- if(ns0>1) rep(1/ns0,ns0)%*%si[s0ind,] else si[s0ind,]
#   ngrad <- dim(btb)[2]
   z <- .Fortran("nlrdtirp",
                 as.double(si),
                 as.integer(ngrad),
                 as.integer(nvox),
                 as.double(btb),
                 as.double(sdcoef),
                 as.double(s0),
                 as.integer(200),
                 as.double(1e-6),
                 res=double((8+ngrad)*nvox),
                 as.integer(8+ngrad),
                 double(ngrad),
                 DUP=FALSE,
                 PACKAGE="dti")$res
    z
}
pnltens <- function(si,grad,s0ind,sdcoef){
#
#  to be used with pmatrix
#
ngrad <- length(grad)/3
lindD <- length(si)/(ngrad+length(s0ind))
dim(si) <- c(ngrad+length(s0ind),lindD)
zmat <- matrix(0,ngrad+8,lindD)
for(i in 1:lindD){
s0 <- si[s0ind,i]
sb <- si[-s0ind,i]
zz <- optim(c(1,0,0,1,0,1),opttensR,method="BFGS",si=sb,s0=s0,grad=grad,sdcoef=sdcoef)
zmat[7,i] <- s0
if(zz$convergence>1&&!is.na(zz$value)&&zz$value<1e12){
zmat[1:6,i] <- rho2D(zz$par)
zmat[8,i] <- zz$value
zmat[-(1:8),i] <- tensRres(zz$par,sb,s0,grad)
} else {
zmat[1:6,i] <- c(1,0,0,1,0,1)
zmat[8,i] <- 1e12
zmat[-(1:8),i] <- rep(0,ngrad)
}
}
zmat
}

pmixtens <- function(x,ngrad0,maxcomp,maxit,pen,grad,reltol,th,penIC,vert){
nvox <- length(x)/(ngrad0+3+maxcomp)
dim(x) <- c(ngrad0+3+maxcomp,nvox)
z <- .C("mixture", 
          as.integer(nvox),
          as.integer(x[(ngrad0+2):(ngrad0+3+maxcomp),]),  
          as.integer(ngrad0),
          as.integer(maxcomp),
          as.integer(maxit),
          as.double(pen),
          as.double(t(grad)),
          as.double(reltol),
          as.double(th),
          as.double(penIC),
          as.double(x[ngrad0+1,]),
          as.double(vert),
          as.double(x[1:ngrad0,]),
          sigma2  = double(nvox),# error variance 
          orient  = double(2*maxcomp*nvox), # phi/theta for all mixture tensors
          order   = integer(nvox),   # selected order of mixture
          lev     = double(2*nvox),         # logarithmic eigenvalues
          mix     = double(maxcomp*nvox),   # mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","lev","mix")]
rbind(z$order,z$sigma2,matrix(z$lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}
pmixtns0 <- function(x,ngrad0,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
nvox <- length(x)/(ngrad0+3+maxcomp)
dim(x) <- c(ngrad0+3+maxcomp,nvox)
z <- .C("mixtrl0", 
          as.integer(nvox),#n1
          as.integer(x[(ngrad0+2):(ngrad0+1+maxcomp),]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad0+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad0,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          mix     = double(maxcomp*nvox),#mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","mix")]
lev <- c(alpha,1)*lambda
rbind(z$order,z$sigma2,matrix(lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}

pmixtns1 <- function(x,ngrad0,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
nvox <- length(x)/(ngrad0+3+maxcomp)
dim(x) <- c(ngrad0+3+maxcomp,nvox)
z <- .C("mixtrl1", 
          as.integer(nvox),#n1
          as.integer(x[(ngrad0+2):(ngrad0+1+maxcomp),]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad0+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad0,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          lambda  = double(nvox),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvox),#mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","lambda","mix")]
lev <- c(alpha,1)*z$lambda
rbind(z$order,z$sigma2,matrix(lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}

pmixtns2 <- function(x,ngrad0,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
nvox <- length(x)/(ngrad0+3+maxcomp)
dim(x) <- c(ngrad0+3+maxcomp,nvox)
z <- .C("mixtrl2", 
          as.integer(nvox),#n1
          as.integer(x[(ngrad0+2):(ngrad0+1+maxcomp),]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad0+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad0,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          alpha   = double(nvox),#alpha_ret alpha=(lambda_1-lambda_2)/lambda_2 
          lambda  = double(nvox),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvox),#mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")] 
lev <- c(z$alpha,1)*z$lambda
rbind(z$order,z$sigma2,matrix(lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}
pgetsii30 <- function(x,maxcomp,dgrad,th,isample0,nsi,nth,nvico,nguess){
         nvox <- length(x)/(nsi+2)
         dim(x) <- c(nsi+2,nvox)
z <- .Fortran("pgtsii30",
         as.double(x[1:nsi,]),
         as.double(x[nsi+1,]),#sigma2
         as.integer(nsi),
         as.integer(nvox),
         as.integer(maxcomp),
         as.double(dgrad),
         as.integer(nvico),
         as.double(th),
         as.integer(nth),
         as.integer(x[nsi+2,]),#indth
         double(nsi*nvico),
         as.integer(isample0),
         as.integer(nguess),
         double(nsi),
         double(nsi*(maxcomp+2)),
         siind=integer((maxcomp+2)*nvox),
         krit=double(nvox),
         as.integer(maxcomp+2),
         DUP=FALSE,
         PACKAGE="dti")[c("siind","krit")]
         dim(z$siind) <- c(maxcomp+2,nvox)
         rbind(z$krit,z$siind)
}
pgetsii31 <- function(x,maxcomp,dgrad,th,isample1,nsi,nth,nvico,nguess,dgradi,maxc){
         nvox <- length(x)/(nsi+3)
         dim(x) <- c(nsi+3,nvox)
z <- .Fortran("pgtsii31",
         as.double(x[1:nsi,]),
         as.double(x[nsi+1,]),#sigma2
         as.integer(nsi),
         as.integer(nvox),
         as.integer(maxcomp),
         as.double(dgrad),
         as.integer(nvico),
         as.integer(x[nsi+3,]),#iandir
         as.double(th),
         as.integer(nth),
         as.integer(x[nsi+2,]),#indth
         double(nsi*nvico),
         as.integer(isample1),
         as.integer(nguess),
         double(nsi),
         double(nsi*(maxcomp+2)),
         siind=integer((maxcomp+2)*nvox),
         krit=double(nvox),
         as.integer(maxcomp+2),
         as.double(dgradi),
         as.double(maxc),
         DUP=FALSE,
         PACKAGE="dti")[c("siind","krit")]
         dim(z$siind) <- c(maxcomp+2,nvox)
         rbind(z$krit,z$siind)
}