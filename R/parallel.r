set.mc.cores <- function(mc.cores=2){
   mc.cores <- min(mc.cores,detectCores())
   num_threads <- as.integer(Sys.getenv("OMP_NUM_THREADS", unset=0))
   if(num_threads <= 0) Sys.setenv("OMP_NUM_THREADS"=mc.cores)
   options("mc.cores"=mc.cores)
   invisible(mc.cores)
}

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
            mccollect(children(jobs), FALSE)
            mckill(children(jobs), if (is.integer(mc.cleanup)) 
                mc.cleanup
            else 15L)
            mccollect(children(jobs))
        }
        if (length(jobs)) {
            mccollect(children(jobs), FALSE)
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
pmatrix <- function (v, FUN, ..., mc.set.seed = TRUE, mc.silent = FALSE, 
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
    l <- if (n <= cores) 
        as.list(data.frame(v))
    else {
        il <- as.integer(n/cores)
        xc <- n - il * cores
        sl <- rep(il, cores)
        if (xc) 
            sl[1:xc] <- il + 1L
        si <- cumsum(c(1L, sl))
        se <- si + c(sl, 0L) - 1L
        lapply(seq_len(cores), function(ix) v[,si[ix]:se[ix]])
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
    if (length(res)<n || length(res)%%n != 0) 
        warning("some results may be missing, folded or caused an error")
    res
}


pdti3Dev <- function(D){
nvox <- dim(D)[2]
.Fortran("dti3Dev",
         as.double(D),
         as.integer(nvox),
         as.logical(rep(TRUE,nvox)),
         ev=double(3*nvox),
         DUP=FALSE,
         PACKAGE="dti")$ev
}

pdtiind3D <- function(D){
nvox <- dim(D)[2]
.Fortran("dtiind3p",
         as.double(D),
         as.integer(nvox),
         res=double(9*nvox),
         DUP=FALSE,
         PACKAGE="dti")$res
}


pnlrdtirg <- function(si,btb,sdcoef,s0ind){
   ns0 <- length(s0ind)
   s0 <- if(ns0>1) rep(1/ns0,ns0)%*%si[s0ind,] else si[s0ind,]
   nvox <- length(s0)
   ngrad <- dim(btb)[2]
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
}

