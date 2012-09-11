pdataframe <- function (v, FUN, ..., mc.set.seed = TRUE, mc.silent = FALSE, 
    mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE) 
{
    require(parallel)
    if (!is.data.frame(v)) 
        stop("'v' must be a dataframe")
    env <- parent.frame()
    cores <- as.integer(mc.cores)
    if (cores < 1L) 
        stop("'mc.cores' must be >= 1")
    if (cores == 1L) 
        return(FUN(v, ...))
    if (mc.set.seed) 
        mc.reset.stream()
    n <- length(v)
    l <- if (n <= cores) 
        as.list(v)
    else {
        il <- as.integer(n/cores)
        xc <- n - il * cores
        sl <- rep(il, cores)
        if (xc) 
            sl[1:xc] <- il + 1L
        si <- cumsum(c(1L, sl))
        se <- si + c(sl, 0L) - 1L
        lapply(seq_len(cores), function(ix) v[si[ix]:se[ix]])
    }
    jobs <- NULL
    cleanup <- function() {
        if (length(jobs) && mc.cleanup) {
            mccollect(parallel:::children(jobs), FALSE)
            parallel:::mckill(parallel:::children(jobs), if (is.integer(mc.cleanup)) 
                mc.cleanup
            else 15L)
            mccollect(parallel:::children(jobs))
        }
        if (length(jobs)) {
            mccollect(parallel:::children(jobs), FALSE)
        }
    }
    on.exit(cleanup())
    FUN <- match.fun(FUN)
    jobs <- lapply(seq_len(min(n, cores)), function(i) mcparallel(FUN(l[[i]], 
        ...), name = i, mc.set.seed = mc.set.seed, silent = mc.silent))
    res <- mccollect(jobs)
    names(res) <- NULL
    res <- do.call(c, res)
    if (length(res) != n) 
        warning("some results may be missing, folded or caused an error")
    res
}
pmatrix <- function(x, FUN, ..., mc.cores = getOption("mc.cores", 2L)){
     require(parallel)
     cl <- makeCluster(mc <- mc.cores)
#
#   parCapply does not pass dimension attributes to FUN !!!!
#
     z <- parCapply(cl, x, FUN, ...)
     stopCluster(cl)
     z
}
plmatrix <- function(x, FUN, ..., mc.cores = getOption("mc.cores", 2L)){
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
pmatrix0 <- function (v, FUN, ..., mc.set.seed = TRUE, mc.silent = FALSE, 
    mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE) 
{
    require(parallel)
    if (!is.matrix(v)) 
        stop("'v' must be a matrix")
    env <- parent.frame()
    cores <- as.integer(mc.cores)
    if (cores < 1L) 
        stop("'mc.cores' must be >= 1")
    if (cores == 1L) 
        return(FUN(v, ...))
    if (mc.set.seed) 
        mc.reset.stream()
    n <- dim(v)[2]
    if (n <= cores)  return(FUN(v, ...))
    il <- as.integer(n/cores)
    xc <- n - il * cores
    sl <- rep(il, cores)
    if (xc) sl[1:xc] <- il + 1L
    si <- cumsum(c(1L, sl))
    se <- si + c(sl, 0L) - 1L
    l <- lapply(seq_len(cores), function(ix) v[,si[ix]:se[ix],drop=FALSE])
    jobs <- NULL
    cleanup <- function() {
        if (length(jobs) && mc.cleanup) {
            mccollect(parallel:::children(jobs), FALSE)
            parallel:::mckill(parallel:::children(jobs), if (is.integer(mc.cleanup)) 
                mc.cleanup
            else 15L)
            mccollect(parallel:::children(jobs))
        }
        if (length(jobs)) {
            mccollect(parallel:::children(jobs), FALSE)
        }
    }
    on.exit(cleanup())
    FUN <- match.fun(FUN)
    jobs <- lapply(seq_len(min(n, cores)), function(i) mcparallel(FUN(l[[i]], 
        ...), name = i, mc.set.seed = mc.set.seed, silent = mc.silent))
    res <- mccollect(jobs)
    names(res) <- NULL
    res <- do.call(c, res)
    if (length(res)<n || length(res)%%n != 0) 
        warning("some results may be missing, folded or caused an error")
    res
}

pdti3Dev <- function(D){
nvox <- length(D)/6
.Fortran("dti3Devp",
         as.double(D),
         as.integer(nvox),
         ev=double(3*nvox),
         DUP=FALSE,
         PACKAGE="dti")$ev
}

pdtiind3D <- function(D){
nvox <- length(D)/6
.Fortran("dtiind3p",
         as.double(D),
         as.integer(nvox),
         res=double(9*nvox),
         DUP=FALSE,
         PACKAGE="dti")$res
}

pdti3Dand <- function(D){
nvox <- length(D)/6
.Fortran("dti3Danp",
         as.double(D),
         as.integer(nvox),
         andir=double(9*nvox),
         DUP=FALSE,
         PACKAGE="dti")$andir
}

pdti3Dall <- function(D){
nvox <- length(D)/6
.Fortran("dti3Dalp",
         as.double(D),
         as.integer(nvox),
         ergs=double(9*nvox),
         DUP=FALSE,
         PACKAGE="dti")$ergs
}

pnlrdtirg <- function(si,btb,sdcoef,s0ind,ngrad){
   ns0 <- length(s0ind)
   nvox <- length(si)/ngrad
   dim(si) <- c(ngrad,length(si)/ngrad)
   s0 <- if(ns0>1) rep(1/ns0,ns0)%*%si[s0ind,] else si[s0ind,]
#   ngrad <- dim(btb)[2]
   z <- .Fortran("nlrdtirp",
                 as.integer(si),
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

pmixtens <- function(x,meth,optmeth,ngrad0,maxcomp,maxit,pen,grad,reltol,th,penIC,vert){
nvox <- length(x)/(ngrad0+3+maxcomp)
dim(x) <- c(ngrad0+3+maxcomp,nvox)
z <- .C("mixture2", 
          as.integer(meth),
          as.integer(optmeth), 
          as.integer(1), # n1
          as.integer(1), # n2
          as.integer(nvox),
          as.integer(rep(1L,nvox)), 
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
#          as.double(orient),
          as.double(t(x[1:ngrad0,])),
          sigma2  = double(nvox),# error variance 
          orient  = double(2*maxcomp*nvox), # phi/theta for all mixture tensors
          order   = integer(nvox),   # selected order of mixture
          lev     = double(2*nvox),         # logarithmic eigenvalues
          mix     = double(maxcomp*nvox),   # mixture weights
          DUPL=FALSE, PACKAGE="dti")[c("sigma2","orient","order","lev","mix")]
          rbind(z$order,z$sigma2,matrix(z$lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
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