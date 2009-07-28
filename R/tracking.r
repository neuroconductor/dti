################################################################
#                                                              #
# Section for tracking() functions (public)                    #
#                                                              #
################################################################

tracking <- function(obj,  ...) cat("Fiber tracking not implemented for this class:",class(obj),"\n")

setGeneric("tracking", function(obj,  ...) standardGeneric("tracking"))

setMethod("tracking","dtiTensor", function(obj, xind=NULL, yind=NULL, zind=NULL, method="LINEPROP", minanindex=0.3, maxangle=30, bgcolor="black", size=1, title=FALSE, ...)
{

  method <- switch(method, "LINEPROP" = 1,
                           1)
  dimx <- obj@ddim[1]
  dimy <- obj@ddim[2]
  dimz <- obj@ddim[3]
  if (is.null(xind)) xind <- 1:dimx
  if (is.null(yind)) yind <- 1:dimy
  if (is.null(zind)) zind <- 1:dimz
  roixa <- min(xind); # this is probably not sufficient
  roixe <- max(xind); # this is probably not sufficient
  roiya <- min(yind); # this is probably not sufficient
  roiye <- max(yind); # this is probably not sufficient
  roiza <- min(zind); # this is probably not sufficient
  roize <- max(zind); # this is probably not sufficient

  dtind <- dtiIndices(obj);

  dd <- .Call("interface_tracking",
              as.double(dtind@andir),
              as.double(dtind@fa),
              as.integer(dimx),
              as.integer(dimy),
              as.integer(dimz),
              as.integer(roixa),
              as.integer(roixe),
              as.integer(roiya),
              as.integer(roiye),
              as.integer(roiza),
              as.integer(roize),
              as.double(obj@voxelext[1]),
              as.double(obj@voxelext[2]),
              as.double(obj@voxelext[3]),
#             as.double(minanindex), # not yet used
#             as.double(maxangle),   # not yet used
#             as.integer(method),    # not yet used (for tracking method)
              DUP=FALSE)

  dim(dd) <- c(length(dd)/6,6);
  
  require(rgl);
  open3d();
  rgl.bg(color=bgcolor)
  rgl.lines(dd[,1],dd[,2],dd[,3],
            color=rgb(abs(dd[,4]),abs(dd[,5]),abs(dd[,6])),
            size=size)
  invisible(dd)
})


