setClass("dti",
         representation(btb    = "matrix",
                        ngrad  = "integer", # = dim(btb)[2]
                        ddim   = "integer",
                        ddim0  = "integer",
                        xind   = "integer",
                        yind   = "integer",
                        zind   = "integer",
                        voxelext = "numeric",
                        source = "character")
         )

setClass("dtiData",
         representation(level = "numeric"),
         contains=c("list","dti"),
         validity=function(object){
          if (sum(c("s0","si") %in% names(object)) != 2) {
            cat("s0 and/or si not in data\n")
            return(invisible(FALSE))
          }
          if (!("array" %in% class(object$s0))) {
            cat("s0 not an array\n")
            return(invisible(FALSE))
          }
          if (!("array" %in% class(object$si))) {
            cat("s0 not an array\n")
            return(invisible(FALSE))
          }
          if (length(dim(object$s0)) != 3) {
            cat("dimension of s0 is",dim(object$s0),", but we want 3 dimensions\n")
            return(invisible(FALSE))
          }
          if (length(dim(object$si)) != 4) {
            cat("dimension of si is",dim(object$si),", but we want 4 dimensions\n")
            return(invisible(FALSE))
          }
          if (any(dim(object$s0) != object@ddim)) {
            cat("dimension of s0 is",dim(object$s0),", but we want",object@ddim,"\n")
            return(invisible(FALSE))
          }
          if (any(dim(object$si) != c(object@ddim,object@ngrad))) {
            cat("dimension of si is",dim(object$si),", but we want",c(object@ddim,object@ngrad),"\n")
            return(invisible(FALSE))
          }
         }
         )
setClass("dtiTensor",
         contains=c("list","dti"),
         validity=function(object){
          if (sum(c("theta","sigma","scorr") %in% names(object)) != 3) {
            cat("s0 and/or si not in data\n")
            return(invisible(FALSE))
          }
          if (!("array" %in% class(object$theta))) {
            cat("s0 not an array\n")
            return(invisible(FALSE))
          }
          if (!("array" %in% class(object$sigma))) {
            cat("s0 not an array\n")
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
          if (sum(c("fa","ra","trc","lambda","eigenv") %in% names(object)) != 5) {
            cat("fa,ra, trc, lambda, or eigenv not in data\n")
            return(invisible(FALSE))
          }
         }
        )
