mfunpl1 <- function(par,siq,grad,pen=1e2){
#
#   evaluate rss for Mixtensor-model (with isotropic component)
#
lpar <- length(par)
m <- (lpar+1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl1",as.double(par),#par(lpar)
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
gmfunpl1 <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar+1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl1g",as.double(par),#par
                as.double(siq),#siq
                as.double(t(grad)),#grad
                as.integer(m),#m
                as.integer(lpar),#lpar
                as.integer(ngrad),#ngrad
                double(ngrad*m),#z
                double(m),#w
                double(ngrad),#work
                double(1),#erg
                double(ngrad),#fv
                double(lpar*ngrad),#dfv
                double(m*ngrad),#qv
                double(2*m*ngrad),#dqv
                gradient=double(lpar),#dh
                PACKAGE="dti")$gradient
}
mfunplwghts1 <- function(par,siq,grad,pen=1e2){
#
#   get weights for Mixtensor-model (with isotropic component) and extract parameters
#
lpar <- length(par)
m <- (lpar+1)/2
ngrad <- dim(grad)[1]
w<-.Fortran("mfunpl1",as.double(par),#par(lpar)
                as.double(siq),#siq(ngrad)
                as.double(t(grad)),#grad(3,ngrad)
                as.integer(m),#number of components
                as.integer(lpar),#number of parameters
                as.integer(ngrad),#number of gradients
                as.double(pen),#penalty for negative weights
                double(ngrad*m),#z(ngrad,m) working array
                w=double(ngrad),#w(ngrad) working array
                double(1),#residual sum of squares
                PACKAGE="dti")$w[1:m]
           wm1 <- w[-1] # first comonent corresponds to isotropic part
           o <- order(wm1,decreasing=TRUE)
           ord <- sum(wm1>0)
           if(ord<m-1){
              o <- o[1:ord]
           }
           sw <- sum(w[w>0])
           lev <- c(par[1],-log(sw))
           if(ord>0){
           mix <- wm1[o]/sw
           } else {
           mix <- NULL
           } 
           or <- matrix(par[2:lpar],2,m-1)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par)
}
mfunpl <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model
#
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
gmfunpl <- function(par,siq,grad){
#
#   evaluate rss for Mixtensor-model
#
lpar <- length(par)
m <- (lpar-1)/2
ngrad <- dim(grad)[1]
.Fortran("mfunplgr",as.double(par),#par
                as.double(siq),#siq
                as.double(t(grad)),#grad
                as.integer(m),#m
                as.integer(lpar),#lpar
                as.integer(ngrad),#ngrad
                double(ngrad*m),#z
                double(m),#w
                double(ngrad),#work
                double(1),#erg
                double(ngrad),#fv
                double(lpar*ngrad),#dfv
                double(m*ngrad),#qv
                double(2*m*ngrad),#dqv
                gradient=double(lpar),#dh
                PACKAGE="dti")$gradient
}
mfunplwghts <- function(par,siq,grad){
#
#   get weights for Mixtensor-model and extract parameters
#
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
           lev <- c(par[1],-log(sw))
           mix <- (w[o]/sw)
           or <- matrix(par[2:lpar],2,m)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par)
}
mfunpl2 <- function(par,siq,grad,pex){
#
#   evaluate rss for Jian-model
#
lpar <- length(par)
m <- (lpar-2)/2
ngrad <- dim(grad)[1]
.Fortran("mfunpl2",as.double(par),
                as.double(siq),
                as.double(t(grad)),
                as.integer(m),
                as.double(pex),
                as.integer(lpar),
                as.integer(ngrad),
                double(ngrad*m),
                double(m),
                double(ngrad),
                erg = double(1),
                PACKAGE="dti")$erg
}
mfunpl2wghts <- function(par,siq,grad,pex){
#
#   get weights for Jian-model and extract parameters
#
lpar <- length(par)
m <- (lpar-2)/2
ngrad <- dim(grad)[1]
w<-.Fortran("mfunpl2",as.double(par),
                as.double(siq),
                as.double(t(grad)),
                as.integer(m),
                as.double(pex),
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
           lev <- exp(c(par[1],par[2]))
           mix <- if(ord==0) NULL else w[o]
           or <- matrix(par[3:lpar],2,m)[,o,drop=FALSE]
           or[1,or[1,]<0] <- or[1,or[1,]<0]+pi
           or[1,or[1,]>pi] <- or[1,or[1,]>pi]-pi
           or[2,or[2,]<0] <- or[2,or[2,]<0]+2*pi
           or[2,or[2,]>2*pi] <- or[2,or[2,]>2*pi]-2*pi
           par <- c(par[1],par[2],or[,1:ord])
list(ord=ord,lev=lev,mix=mix,orient=or,par=par)
}


dwiMixtensor.new <- function(object, ...) cat("No dwiMixtensor.new calculation defined for this class:",class(object),"\n")

setGeneric("dwiMixtensor.new", function(object,  ...) standardGeneric("dwiMixtensor.new"))

setMethod("dwiMixtensor.new",
          "dtiData",
          function(object,
                   maxcomp = 3,
                   maxneighb = 7,
                   reltol = 1e-8,
                   maxit = 5000,
                   optmethod = "Nelder-Mead")
           {
             cat("entering method dwiMixtensor\n")

             if (maxcomp > 4) stop("maximum number of tensor components is 4\n") # otherwise change number of parameters in C code!

             args <- sys.call(-1)
             args <- c(object@call, args)
             ngrad <- object@ngrad
             ddim <- object@ddim
             s0ind <- object@s0ind
             ns0 <- length(s0ind)

             cat("determine outliers ... ")
             z <- .Fortran("outlier",
                           as.double(object@si),
                           as.integer(prod(ddim)),
                           as.integer(ngrad),
                           as.logical((1:ngrad)%in%s0ind),
                           as.integer(ns0),
                           si      = integer(prod(ddim)*ngrad),
                           index   = integer(prod(ddim)),
                           lindex  = integer(1),
                           DUP     = FALSE,
                           PACKAGE = "dti")[c("si","index","lindex")]
             si <- array(z$si, c(ddim, ngrad))
             index <- if (z$lindex>0) z$index[1:z$lindex] else numeric(0)
             rm(z)
             cat("done\n")

             cat("prepare data and initial estimates ... ")
             # prepare data for optim
             s0 <- si[,,,s0ind,drop=FALSE]
             if (length(s0ind)>1) s0 <- apply(s0, 1:3, mean) else dim(s0) <- dim(s0)[1:3]
             # normalized DW data
             mask <- s0 > object@level
             siq <- si[,,,-s0ind,drop=FALSE]
             dim(siq) <- c(prod(ddim),ngrad-ns0)
             siq[mask,] <- sweep(siq[mask,],1,s0[mask],"/")
             dim(siq) <- c(ddim,ngrad-ns0)
             # heuristics to avoid DWI that are larger than s0.
             siqmed <- apply(siq, 1:3, median)
             siqmed[siqmed < .9] <- .9
             siqmed[siqmed > .99] <- .99
             siq <- sweep(siq,1:3,siqmed,pmin)
             # mask for calculation
             mask <- s0 > object@level
             grad <- t(object@gradient[,-s0ind])
             siind <- getsiind(siq, mask, grad, maxcomp, maxc=.866)
             # siind[1,,,] contains number of potential directions 
             # siind[-1,,,] contains indices of grad corresponding to these directions
             cat("done\n")

             cat("optimizing ... ")
             mm <- switch(optmethod, "Nelder-Mead" = 1,
                                     "BFGS"        = 3,
                                     1)
             pl <- 1

             # perform voxelwise optimization and order selection of tensor mixture model
             a <- .C("mixturepl",
                     as.integer(mm),                         # select mixture method
                     as.integer(prod(ddim)),                 # number of voxels
                     as.integer(mask),                       # calculation mask
                     as.double(siq),                         # DWI without s0
                     as.integer(siind),                      # DWI indices of local minima
                     as.integer(ngrad - length(s0ind)),      # number of DWI
                     as.double(grad),                        # gradient directions
                     as.integer(maxcomp),                    # max number of gradient neighbors
                     as.integer(pl),                         # exp for Jian model
                     as.integer(maxit),                      # max number of iterations for optim
                     as.double(reltol),                      # reltol crit for optim
                     order   = integer(prod(ddim)),           # selected order of mixture
                     lev     = double(2*prod(ddim)),         # logarithmic eigenvalues
                     mix     = double(maxcomp*prod(ddim)),   # mixture weights
                     orient  = double(2*maxcomp*prod(ddim)), # phi/theta for all mixture tensors
                     sigma2  = double(prod(ddim)),           # parameter variance ???
                     DUP     = FALSE,
                     PACKAGE = "dti")[c("order", "lev", "mix", "orient", "sigma2")]

             # set dimension attr
             dim(a$order) <- ddim;
             dim(a$lev) <- c(2, ddim);
             dim(a$mix) <- c(maxcomp, ddim);
             dim(a$orient) <- c(2, maxcomp, ddim);
             dim(a$sigma2) <- ddim;
             cat("done\n")

             # create and return new object
             invisible(new("dwiMixtensor",
                           call        = args,
                           ev          = a$lev,
                           mix         = a$mix,
                           orient      = a$orient,
                           order       = a$order,
                           p           = pl,
                           th0         = s0,
                           sigma       = a$sigma2,
                           scorr       = array(1, c(1,1,1)), # ???
                           bw          = c(0,0,0),           # ???
                           mask        = mask,
                           hmax        = 1,                  # ???
                           gradient    = object@gradient,
                           btb         = object@btb,
                           ngrad       = object@ngrad,
                           s0ind       = object@s0ind,
                           replind     = object@replind,
                           ddim        = object@ddim,
                           ddim0       = object@ddim0,
                           xind        = object@xind,
                           yind        = object@yind,
                           zind        = object@zind,
                           voxelext    = object@voxelext,
                           level       = object@level,
                           orientation = object@orientation,
                           source      = object@source,
                           outlier     = index,
                           scale       = 1,
                           method      = "mixtensor")
                       )
           })
