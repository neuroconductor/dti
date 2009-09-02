################################################################
#                                                              #
# Section for tracking() functions (public)                    #
#                                                              #
################################################################

tracking <- function(obj,  ...) cat("Fiber tracking not implemented for this class:",class(obj),"\n")

setGeneric("tracking", function(obj,  ...) standardGeneric("tracking"))

setMethod("tracking","dtiTensor", function(obj, xind=NULL, yind=NULL, zind=NULL, method="LINEPROP", minanindex=0.3, maxangle=30)
{

  args <- sys.call(-1)
  args <- c(obj@call,args)
  imethod <- switch(method, "LINEPROP" = 1,
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
  if(sum(dtind@fa[xind,yind,zind]>minanindex)==0){
     cat("No fiber with sufficint FA in region of interest\n")
     return(invisible(FALSE))
  }

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
                roix   = as.integer(range(xind)),
                roiy   = as.integer(range(yind)),
                roiz   = as.integer(range(zind)),
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

setMethod("tracking","dtiIndices", function(obj, xind=NULL, yind=NULL, zind=NULL, method="LINEPROP", minanindex=0.3, maxangle=30)
{

  args <- sys.call(-1)
  args <- c(obj@call,args)
  imethod <- switch(method, "LINEPROP" = 1,
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
  
  if(sum(obj@fa[xind,yind,zind]>minanindex)==0){
     cat("No fiber with sufficint FA in region of interest\n")
     return(invisible(FALSE))
  }

  dd <- .Call("interface_tracking",
              as.double(obj@andir),
              as.double(obj@fa),
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
                roix   = as.integer(range(xind)),
                roiy   = as.integer(range(yind)),
                roiz   = as.integer(range(zind)),
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
