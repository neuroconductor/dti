################################################################
#                                                              #
# Section for tracking() functions (public)                    #
#                                                              #
################################################################

tracking <- function(obj,  ...) cat("Fiber tracking not implemented for this class:",class(obj),"\n")

setGeneric("tracking", function(obj,  ...) standardGeneric("tracking"))

setMethod("tracking","dtiTensor", function(obj, roix=NULL, roiy=NULL, roiz=NULL, method="LINEPROP", minanindex=0.3, maxangle=30, subsample=1)
{

  args <- sys.call(-1)
  args <- c(obj@call,args)
  imethod <- switch(method, "LINEPROP" = 1,
                           1)
  dimx <- obj@ddim[1]
  dimy <- obj@ddim[2]
  dimz <- obj@ddim[3]
  if (is.null(roix)) roix <- 1:dimx
  if (is.null(roiy)) roiy <- 1:dimy
  if (is.null(roiz)) roiz <- 1:dimz
  roixa <- min(roix); # this is probably not sufficient
  roixe <- max(roix); # this is probably not sufficient
  roiya <- min(roiy); # this is probably not sufficient
  roiye <- max(roiy); # this is probably not sufficient
  roiza <- min(roiz); # this is probably not sufficient
  roize <- max(roiz); # this is probably not sufficient

  dtind <- dtiIndices(obj);
  if(sum(dtind@fa[roix,roiy,roiz]>minanindex)==0){
     cat("No fiber with sufficint FA in region of interest\n")
     return(invisible(FALSE))
  }

  andir <- dtind@andir
  fa <- dtind@fa

  if ((subsample != as.integer(subsample)) | (subsample < 1)) subsample <- 1
  if (subsample > 1) {
    indx <- rep(1:dimx, rep(subsample,dimx))
    indy <- rep(1:dimy, rep(subsample,dimy))
    indz <- rep(1:dimz, rep(subsample,dimz))
    fa <- fa[indx, indy, indz]
    andir <- andir[,indx, indy, indz]
    dimx <- subsample*dimx
    dimy <- subsample*dimy
    dimz <- subsample*dimz
    roixa <- (roixa-1)*subsample+1
    roixe <- roixe*subsample
    roiya <- (roiya-1)*subsample+1
    roiye <- roiye*subsample
    roiza <- (roiza-1)*subsample+1
    roize <- roize*subsample
  }
  
  
  dd <- .Call("interface_tracking",
              as.double(andir),
              as.double(fa),
              as.integer(dimx),
              as.integer(dimy),
              as.integer(dimz),
              as.integer(roixa),
              as.integer(roixe),
              as.integer(roiya),
              as.integer(roiye),
              as.integer(roiza),
              as.integer(roize),
              as.double(obj@voxelext[1]/subsample),
              as.double(obj@voxelext[2]/subsample),
              as.double(obj@voxelext[3]/subsample),
              as.double(minanindex), # not yet used
              as.double(maxangle),   # not yet used
#             as.integer(imethod),    # not yet used (for tracking method)
              DUP=FALSE)

  dim(dd) <- c(length(dd)/6,6);
  dd <- reduce.fibers(dd)
  istartfiber <- ident.fibers(dd)
  invisible(new("dwiFiber",
                call  = args,
                fibers = dd,
                startind = as.integer(istartfiber),
                roix   = as.integer(range(roix)),
                roiy   = as.integer(range(roiy)),
                roiz   = as.integer(range(roiz)),
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
                source = obj@source,
                method = method,
                minanindex = minanindex,
                maxangle = maxangle)
            )
})

setMethod("tracking","dtiIndices", function(obj, roix=NULL, roiy=NULL, roiz=NULL, method="LINEPROP", minanindex=0.3, maxangle=30, subsample = 1)
{

  args <- sys.call(-1)
  args <- c(obj@call,args)
  imethod <- switch(method, "LINEPROP" = 1,
                           1)

  dimx <- obj@ddim[1]
  dimy <- obj@ddim[2]
  dimz <- obj@ddim[3]
  if (is.null(roix)) roix <- 1:dimx
  if (is.null(roiy)) roiy <- 1:dimy
  if (is.null(roiz)) roiz <- 1:dimz
  roixa <- min(roix); # this is probably not sufficient
  roixe <- max(roix); # this is probably not sufficient
  roiya <- min(roiy); # this is probably not sufficient
  roiye <- max(roiy); # this is probably not sufficient
  roiza <- min(roiz); # this is probably not sufficient
  roize <- max(roiz); # this is probably not sufficient

  if(sum(obj@fa[roix,roiy,roiz]>minanindex)==0){
     cat("No fiber with sufficint FA in region of interest\n")
     return(invisible(FALSE))
  }

  andir <- obj@andir
  fa <- obj@fa

  if ((subsample != as.integer(subsample)) | (subsample < 1)) subsample <- 1
  if (subsample > 1) {
    indx <- rep(1:dimx, rep(subsample,dimx))
    indy <- rep(1:dimy, rep(subsample,dimy))
    indz <- rep(1:dimz, rep(subsample,dimz))
    fa <- fa[indx, indy, indz]
    andir <- andir[,indx, indy, indz]
    dimx <- subsample*dimx
    dimy <- subsample*dimy
    dimz <- subsample*dimz
    roixa <- (roixa-1)*subsample+1
    roixe <- roixe*subsample
    roiya <- (roiya-1)*subsample+1
    roiye <- roiye*subsample
    roiza <- (roiza-1)*subsample+1
    roize <- roize*subsample
  }
  
  
  dd <- .Call("interface_tracking",
              as.double(andir),
              as.double(fa),
              as.integer(dimx),
              as.integer(dimy),
              as.integer(dimz),
              as.integer(roixa),
              as.integer(roixe),
              as.integer(roiya),
              as.integer(roiye),
              as.integer(roiza),
              as.integer(roize),
              as.double(obj@voxelext[1]/subsample),
              as.double(obj@voxelext[2]/subsample),
              as.double(obj@voxelext[3]/subsample),
              as.double(minanindex), # not yet used
              as.double(maxangle),   # not yet used
#             as.integer(imethod),    # not yet used (for tracking method)
              DUP=FALSE)

  dim(dd) <- c(length(dd)/6,6);
#  dd <- reduce.fibers(dd)
  istartfiber <- ident.fibers(dd)
  invisible(new("dwiFiber",
                call  = args,
                fibers = dd,
                startind = as.integer(istartfiber),
                roix   = as.integer(range(roix)),
                roiy   = as.integer(range(roiy)),
                roiz   = as.integer(range(roiz)),
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
                source = obj@source,
                method = method,
                minanindex = minanindex,
                maxangle = maxangle)
            )
})

setMethod("tracking", "dwiMixtensor",
function(obj, roix=NULL, roiy=NULL, roiz=NULL, method="LINEPROP", minanindex=0.3, maxangle=30, subsample = 1)
{

  args <- sys.call(-1)
  args <- c(obj@call,args)
  imethod <- switch(method, "LINEPROP" = 1,
                           1)

  dimx <- obj@ddim[1]
  dimy <- obj@ddim[2]
  dimz <- obj@ddim[3]
  if (is.null(roix)) roix <- 1:dimx
  if (is.null(roiy)) roiy <- 1:dimy
  if (is.null(roiz)) roiz <- 1:dimz
  roixa <- min(roix); # this is probably not sufficient
  roixe <- max(roix); # this is probably not sufficient
  roiya <- min(roiy); # this is probably not sufficient
  roiye <- max(roiy); # this is probably not sufficient
  roiza <- min(roiz); # this is probably not sufficient
  roize <- max(roiz); # this is probably not sufficient

  ex <- extract(obj, c("andir", "order", "gfa", "mix"))

  if(sum(ex@gfa[roix,roiy,roiz] > minanindex)==0){
     cat("No fiber with sufficint FA in region of interest\n")
     return(invisible(FALSE))
  }

  if ((subsample != as.integer(subsample)) | (subsample < 1)) subsample <- 1
  if (subsample > 1) {
    indx <- rep(1:dimx, rep(subsample,dimx))
    indy <- rep(1:dimy, rep(subsample,dimy))
    indz <- rep(1:dimz, rep(subsample,dimz))
    ex$fa <- ex$fa[indx, indy, indz]
    ex$order <- ex$order[indx, indy, indz]
    ex$andir <- ex$andir[,,indx, indy, indz, drop=FALSE]
    dimx <- subsample*dimx
    dimy <- subsample*dimy
    dimz <- subsample*dimz
    roixa <- (roixa-1)*subsample+1
    roixe <- roixe*subsample
    roiya <- (roiya-1)*subsample+1
    roiye <- roiye*subsample
    roiza <- (roiza-1)*subsample+1
    roize <- roize*subsample
  }
  maxorder <- dim(ex$andir)[2]
  
  dd <- .Call("interface_tracking_mixtensor",
              as.double(ex$andir), # dim = c(3, maxorder, dimx, dimy, dimz)
              as.double(ex$order), # NEW! dim = c(dimx, dimy, dimz)
              as.double(ex$gfa),    # dim = c(dimx, dimy, dimz)
              as.double(ex$mix),    # NEW! dim = c(maxorder, dimx, dimy, dimz)
              as.integer(maxorder), # NEW!
              as.integer(dimx),
              as.integer(dimy),
              as.integer(dimz),
              as.integer(roixa),
              as.integer(roixe),
              as.integer(roiya),
              as.integer(roiye),
              as.integer(roiza),
              as.integer(roize),
              as.double(obj@voxelext[1]/subsample),
              as.double(obj@voxelext[2]/subsample),
              as.double(obj@voxelext[3]/subsample),
              as.double(minanindex),
              as.double(maxangle),
#             as.integer(imethod),    # not yet used (for tracking method)
              DUP=FALSE)

  dim(dd) <- c(length(dd)/6,6);
#  dd <- reduce.fibers(dd)
  istartfiber <- ident.fibers(dd)
  invisible(new("dwiFiber",
                call  = args,
                fibers = dd,
                startind = as.integer(istartfiber),
                roix   = as.integer(range(roix)),
                roiy   = as.integer(range(roiy)),
                roiz   = as.integer(range(roiz)),
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
                source = obj@source,
                method = method,
                minanindex = minanindex,
                maxangle = maxangle)
            )
})

ident.fibers <- function(mat){
#
#  Identify indices in mat where a new fiber starts
#
   dd <- dim(mat)
   if(dd[2]!=6){
      warning("Incorrect dimensions for fiber array")
   }
   dd <- dd[1]
   mat <- mat[,1:3]
   dim(mat) <- c(2,dd/2,3)
   dmat <- mat[2,-(dd/2),]-mat[1,-1,]
   fiberends <- (1:(dd/2-1))[apply(dmat^2,1,sum)>0]
   c(0,fiberends)*2+1
}

reduce.fibers <- function(mat){
   dd <- dim(mat)
#
#  clean up fiber description in dd
#  removes instances in dd that are either redundant or would not show up on display
#
   if(dd[2]!=6){
      warning("Incorrect dimensions for fiber array")
   }
   dd <- dd[1]
   dim(mat) <- c(2,dd/2,6)
   dmat <- mat[1,,1:3]-mat[2,,1:3]
   remove <- apply(dmat^2,1,sum)==0
   if(sum(remove)>0) warning(paste("Found ",sum(remove)," instances, where begin and end of a line segment coincide"))
   mat <- mat[,!remove,]
   dim(mat) <- c(2*dim(mat)[2],6)
   mat
}
