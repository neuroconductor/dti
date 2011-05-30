gethseqse3<-
function (kstar, grad, kexp=.1, vext = c(1, 1), dist = "SE3", 
    verbose = FALSE) 
{
    ngrad <- dim(grad)[2]
    h <- kappa <- vr <- vr0 <- matrix(0,ngrad,kstar)
    prt0 <- Sys.time()
    cat("get sequences of bw, kappa up to kstar=", kstar, " ")
    n <- 0
    for(i in 1:ngrad){
       z <- .Fortran("gethse3i",
                     as.integer(i),
                     as.integer(kstar),
                     as.double(grad),
                     double(3*ngrad),
                     double(2*ngrad),
                     as.integer(ngrad),
                     as.double(kexp),
                     as.double(vext),
                     h=double(kstar),
                     kappa=double(kstar),
                     vr=double(kstar),
                     vr0=double(kstar),
                     n=integer(1),
                     DUPL=FALSE,
                     PACKAGE="dti")[c("h","kappa","vr","vr0","n")]
       h[i,] <- z$h
       kappa[i,] <- z$kappa
       vr[i,] <- z$vr
       vr0[i,] <- z$vr0 
       n <- n+z$n
        if (verbose) {
            cat("i",i,"h[i,.]:", signif(h[i,],3), "\n kappa[i,.]:", signif(kappa[i,],3), "time elapsed:", 
                format(difftime(Sys.time(), prt0), digits = 3), 
                "\n")
        }
        else {
            cat(".")
        }
        gc()
    }
        cat("number of positive weights:",n,"\n")
    list(h=h,kappa=kappa,vred=vr,vred0=vr0,n=n)
}
gethseqfullse3 <-
function (kstar, gradstats, kappa0=NULL, kexp=.1, vext = c(1, 1), dist = "SE3", 
    verbose = FALSE) 
{
    ngrad <- dim(gradstats$bg)[2]
    h <- kappa <- vr <- vr0 <- matrix(0,ngrad,kstar)
    if(is.null(kappa0)) kappa0 <- pi/ngrad
    prt0 <- Sys.time()
    cat("get sequences of bw, kappa up to kstar=", kstar, " ")
    n <- 0
    for(i in 1:ngrad){
       z <- .Fortran("ghfse3i",
                     as.integer(i),#i4
                     as.integer(kstar),#kstar
                     as.double(gradstats$k456),
                     as.double(gradstats$bghat),
                     as.double(gradstats$nbg),
                     as.double(gradstats$nbghat),
                     as.integer(ngrad),
                     as.double(kexp),#kexp,
                     as.double(kappa0),#kappa0
                     as.double(vext),#vext
                     h=double(kstar),
                     kappa=double(kstar),#
                     vr=double(kstar),#
                     vr0=double(kstar),#
                     n=integer(1),#
                     DUPL=FALSE,
                     PACKAGE="dti")[c("h","kappa","vr","vr0","n")]
       h[i,] <- z$h
       kappa[i,] <- z$kappa
       vr[i,] <- z$vr
       vr0[i,] <- z$vr0 
       n <- n+z$n
        if (verbose) {
            cat("i",i,"h[i,.]:", signif(h[i,],3), "\n kappa[i,.]:", signif(kappa[i,],3), "time elapsed:", 
                format(difftime(Sys.time(), prt0), digits = 3), 
                "\n")
        }
        else {
            cat(".")
        }
        gc()
    }
        cat("number of positive weights:",n,"\n")
    list(h=h,kappa=kappa,vred=vr,vred0=vr0,n=n)
}


