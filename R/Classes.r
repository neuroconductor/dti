setClass("dti",
         representation(.Data = "list",
                        btb    = "matrix",
                        ngrad  = "integer", # = dim(btb)[2]
                        s0ind  = "integer", # indices of s0 images
                        replind = "integer", # replications in gradient design
                        ddim   = "integer",
                        ddim0  = "integer",
                        xind   = "integer",
                        yind   = "integer",
                        zind   = "integer",
                        voxelext = "numeric",
                        source = "character"),
         )

setClass("dtiData",
         representation(si = "array", level = "numeric"),
         contains=c("list","dti"),
         validity=function(object){
          if (any(dim(object@si)!=c(object@ddim,object@ngrad))) {
            cat("incorrect dimension of image data\n")
            return(invisible(FALSE))
          }
          if (length(object@s0ind)<1) {
            cat("no S_0 images, parameters not identifiable \n")
             return(invisible(FALSE))
         }
         }
         )
setClass("dtiTensor",
         representation(method = "character",
                        D      = "array",
                        th0    = "array",
                        sigma  = "array",
                        scorr  = "array",
                        bw     = "numeric",
                        mask   = "array",
                        hmax   = "numeric"),
         contains=c("list","dti"),
         validity=function(object){
          if (any(dim(object@D)!=c(6,object@ddim))) {
            cat("invalid dimension of tensor array D \n")
            return(invisible(FALSE))
          }
          if (any(dim(object@th0)!=object@ddim)) {
            cat("invalid dimension of array th0\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@sigma)!=object@ddim)) {
            cat("invalid dimension of array sigma\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@mask)!=object@ddim)) {
            cat("invalid dimension of array mask\n")
            return(invisible(FALSE))
          }
          if (!is.logical(object@mask)) {
            cat("invalid type of array mask, should be logical\n")
            return(invisible(FALSE))
          }
          if (length(dim(object@scorr))!=3) {
            cat("invalid dimension of scorr\n")
            return(invisible(FALSE))
          }
          if (length(object@bw)!=3) {
            cat("invalid length of bw\n")
            return(invisible(FALSE))
          }
          if (!(object@method %in% c("linear","nonlinear","unknown"))) {
            cat("method should specify either linear or nonlinear or unknown\n")
            return(invisible(FALSE))
          }
         }
         )

setClass("dtiIndices",
         representation(fa     = "array",
                        ra     = "array",
                        trc    = "array",
                        bary   = "array",
                        lambda = "array",
                        eigenv = "array"),
         contains=c("list","dti"),
          validity=function(object){
          if (any(dim(object@fa)!=object@ddim)) {
            cat("invalid dimension of array fa\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@ra)!=object@ddim)) {
            cat("invalid dimension of array ra\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@trc)!=object@ddim)) {
            cat("invalid dimension of array trc\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@bary)!=c(object@ddim,3))) {
            cat("invalid dimension of array bary\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@lambda)!=c(object@ddim,3))) {
            cat("invalid dimension of array bary\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@eigenv)!=c(object@ddim,9))) {
            cat("invalid dimension of array eigenv\n")
            return(invisible(FALSE))
          }
         }
        )
