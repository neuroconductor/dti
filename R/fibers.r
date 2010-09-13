selectFibers <- function(obj,  ...) cat("Selection of fibers is not implemented for this class:",class(obj),"\n")

setGeneric("selectFibers", function(obj,  ...) standardGeneric("selectFibers"))

setMethod("selectFibers","dwiFiber", function(obj, roix=NULL, roiy=NULL, roiz=NULL, nroimask=NULL, minlength=1)
{
#)
#     extract fiber information and descriptions
#
  args <- sys.call(-1)
  args <- c(obj@call,args)
  fibers <- obj@fibers
  fiberstart <- obj@startind
  fiberlength <- diff(c(fiberstart,dim(fibers)[1]+1))
  if(minlength>1){
#
#   eliminate fibers shorter than minlength
#
  remove <- (1:length(fiberstart))[fiberlength<minlength]
  if(length(remove)>0){
     inda <- fiberstart[remove]
     inde <- c(fiberstart,dim(fibers)[1]+1)[remove+1]-1
     fibers <- fibers[-c(mapply(":",inda,inde),recursive=TRUE),]
     fiberlength <- fiberlength[-remove]
     fiberstart <- c(0,cumsum(fiberlength))[1:length(fiberlength)]+1
  }
  }
  roimask <- as.integer(obj@roimask)
  mroimask <- max(roimask)
  if(mroimask>127){
     warning("Recursive use of regions of interest is limited to 7")
     roix <- NULL
     roimask <- NULL
  }
if(!(is.null(roix)&&is.null(roiy)&&is.null(roiz)&&is.null(nroimask))){
#
#    no region of interest specified otherwise
#
  if(is.null(nroimask)){ 
     nroimask <- array(0,obj@ddim)
     if(is.null(roix)) roix <- 1:obj@ddim[1]
     if(is.null(roiy)) roiy <- 1:obj@ddim[2]
     if(is.null(roiz)) roiz <- 1:obj@ddim[3]
     nroimask[roix,roiy,roiz] <- 1
  }
  lnewmask <- as.integer(log(mroimask,2))+1
  roimask[nroimask>0] <- roimask[nroimask>0]+2^lnewmask
  cat("dimension of roimask",dim(roimask),"\n")
  z <- .Fortran("roifiber",
                as.double(fibers),
                newfibers=double(prod(dim(fibers))),#new fibers
                as.integer(dim(fibers)[1]),
                integer(3*max(fiberlength)),#array for fiber in vcoord
                as.integer(max(fiberlength)),# maximum fiberlength
                start=as.integer(fiberstart),
                as.integer(fiberlength),
                as.integer(length(fiberstart)),#number of fibers
                as.logical(nroimask>0),#roi
                as.integer(obj@ddim[1]),
                as.integer(obj@ddim[2]),
                as.integer(obj@ddim[3]),
                as.double(obj@voxelext),
                sizenf=integer(1),
                nfiber=integer(1),
                DUP=TRUE,
                PACKAGE="dti")[c("newfibers","sizenf","start","nfiber")]
  if(z$sizenf>1) {
     fibers <- array(z$newfibers,dim(fibers))[1:z$sizenf,]
     fiberstart <- z$start[1:z$nfiber]
  }  else {
     warning("No fibers found, return original")
  }
  }
  invisible(new("dwiFiber",
                call  = args,
                fibers = fibers,
                startind = as.integer(fiberstart),
                roimask = as.raw(roimask),
                gradient = obj@gradient,
                btb   = obj@btb,
                ngrad = obj@ngrad, # = dim(btb)[2]
                s0ind = obj@s0ind,
                replind = obj@replind,
                ddim  = obj@ddim,
                ddim0 = obj@ddim0,
                xind  = obj@xind,
                yind  = obj@yind,
                zind  = obj@zind,
                voxelext = obj@voxelext,
                level = obj@level,
                orientation = obj@orientation,
                rotation = obj@rotation,
                source = obj@source,
                method = obj@method,
                minanindex = obj@minanindex,
                maxangle = obj@maxangle)
            )
})

reduceFibers <- function(fiberobj,  ...) cat("Selection of fibers is not implemented for this class:",class(obj),"\n")

setGeneric("reduceFibers", function(fiberobj,  ...) standardGeneric("reduceFibers"))

setMethod("reduceFibers","dwiFiber", function(fiberobj, maxdist=1)
{
   fiberobj <- sort.fibers(fiberobj)
   fibers <- fiberobj@fibers[,1:3]
   nsegm <- dim(fibers)[1]
   startf <- fiberobj@startind
   endf <- c(startf[-1]-1,nsegm)
   nfibers <- length(startf)
   keep <- .Fortran("reducefi",
                    as.double(t(fibers)),
                    as.integer(nsegm),
                    as.integer(startf),
                    as.integer(endf),
                    as.integer(nfibers),
                    keep=logical(nfibers),
                    as.double(maxdist),
                    DUP=FALSE,
                    PACKAGE="dti")$keep
    startf <- startf[keep]
    endf <- endf[keep]
    ind <- rep(startf,endf-startf+1)+sequence(endf-startf+1)-1
    fiberobj@fibers <- fiberobj@fibers[ind,]
    fiberobj@startind <- c(0,cumsum(endf-startf+1))[1:length(startf)]+1
    fiberobj
}
)

#ident.fibers <- function(mat){
#
#  Identify indices in mat where a new fiber starts
#
#   dd <- dim(mat)
#   if(dd[2]!=3){
#      warning("Incorrect dimensions for fiber array")
#   }
#   dd <- dd[1]
#   z <- .Fortran("fibersta",
#                 as.double(mat),
#                 as.integer(dd/2),
#                 fiberstarts=integer(dd/2),#thats more the maximum needed
#                 nfibers=integer(1),
#                 DUP=FALSE,
#                 PACKAGE="dti")[c("fiberstarts","nfibers")]
#   z$fiberstarts[1:z$nfibers]
#}

sort.fibers <- function(fiberobj){
#
#  sort fiber array such that longest fibers come first
#
   fibers <- fiberobj@fibers
   nfs <- dim(fibers)[1]
   starts <- fiberobj@startind
   ends <- c(starts[-1]-1,nfs)
   fiberlength <- diff(c(starts,nfs+1))
   of <- order(fiberlength,decreasing=TRUE)
   ind <-  rep(starts[of],ends[of]-starts[of]+1)+sequence(ends[of]-starts[of]+1)-1
   fiberobj@fibers <- fiberobj@fibers[ind,]
   fiberobj@startind <- as.integer(c(0,cumsum(ends[of]-starts[of]+1))[1:length(starts)]+1)
   fiberobj
} 

  compactFibers <- function(fibers,startind){
  n <- dim(fibers)[1]
  endind <- c(startind[-1]-1,n)
  ind <- 1:n
  ind <- ind[ind%%2==1 | ind%in%endind]
  list(fibers=fibers[ind,],startind=(startind-1)/2+1:length(startind))
  }

  expandFibers <- function(fibers,startind){
  n <- dim(fibers)[1]
  endind <- c(startind[-1]-1,n)
  ind <- rep(2,max(endind))
  ind[startind] <- 1
  ind[endind] <- 1
  ind <-  rep(startind,2*(endind-startind))+rep(sequence((endind-startind)+1)-1,ind)
  startind[-1] <- startind[-1]+cumsum(diff(startind)-2)
  list(fibers=fibers[ind,],startind=startind)
  }
