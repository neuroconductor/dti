gethseqse3 <-
function (kstar, grad, kappa, vext = c(1, 1), dist = "SE3", 
    verbose = FALSE) 
{
    ngrad <- dim(grad)[2]
    h <- vr <- matrix(0,ngrad,kstar)
    if(is.null(kappa)) kappa <- pi/ngrad
    if(length(kappa)<kstar) kappa <- rep(kappa[1],kstar)
    prt0 <- Sys.time()
    cat("get sequences of bw, kappa up to kstar=", kstar, " ")
    n <- 0
    for(i in 1:ngrad){
       z <- .Fortran("gethse3i",
                     as.integer(i),
                     as.integer(kstar),
                     as.double(grad),
                     as.double(kappa[1]),
                     double(3*ngrad),
                     double(2*ngrad),
                     as.integer(ngrad),
                     as.double(vext),
                     h=double(kstar),
                     vr=double(kstar),
                     n=integer(1),
                     DUPL=FALSE,
                     PACKAGE="dti")[c("h","vr","n")]
       h[i,] <- z$h
       vr[i,] <- z$vr 
       n <- n+z$n
        if (verbose) {
            cat("i",i,"h[i,.]:", signif(h[i,],3), "\n")
        }
        else {
            cat(".")
        }
        gc()
    }
        cat("number of positive weights:",n,"time elapsed:", 
                format(difftime(Sys.time(), prt0), digits = 3),"\n")
    list(h=h,kappa=rep(kappa,kstar),vred=vr,n=n)
}
gethseqfullse3 <-
function (kstar, gradstats, kappa=NULL, vext = c(1, 1), dist = "SE3", 
    verbose = FALSE) 
{
    ngrad <- dim(gradstats$bg)[2]
    h  <- vr  <- matrix(0,ngrad,kstar)
    if(is.null(kappa)) kappa <- pi/ngrad
    if(length(kappa)<kstar) kappa <- rep(kappa[1],kstar)
    prt0 <- Sys.time()
    cat("get sequences of bw, kappa up to kstar=", kstar, " ")
    n <- 0
    for(i in 1:ngrad){
       z <- .Fortran("ghfse3i",
                     as.integer(i),#i4
                     as.integer(kstar),#kstar
                     as.double(gradstats$k456),
                     as.double(gradstats$nbg),
                     as.double(gradstats$nbghat),
                     as.integer(ngrad),
                     as.double(kappa),#kappa
                     as.double(vext),#vext
                     h=double(kstar),
                     vr=double(kstar),#
                     n=integer(1),#
                     DUPL=FALSE,
                     PACKAGE="dti")[c("h","vr","n")]
       h[i,] <- z$h
       vr[i,] <- z$vr 
       n <- n+z$n
        if (verbose) {
            cat("i",i,"h[i,.]:", signif(h[i,],3), "\n")
        }
        else {
            cat(".")
        }
        gc()
    }
        cat("number of positive weights:",n,"time elapsed:", 
                format(difftime(Sys.time(), prt0), digits = 3), 
                "\n")
    list(h=h,kappa=kappa,vred=vr,n=n)
}


