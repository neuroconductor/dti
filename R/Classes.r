setClass("dti",
         representation(btb    = "matrix",
                        ngrad  = "integer", # = dim(btb)[2]
                        s0ind  = "integer", # indices of s0 images
                        replind = "integer", # replications in gradient design
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
          if (!("si") %in% names(object)) {
            cat("si not in data\n")
            return(invisible(FALSE))
          }
          if (!("array" %in% class(object$si))) {
            cat("s0 not an array\n")
            return(invisible(FALSE))
          }
          if (length(dim(object$si)) != 4) {
            cat("dimension of si is",dim(object$si),", but we want 4 dimensions\n")
            return(invisible(FALSE))
          }
          if (any(dim(object$si) != c(object@ddim,object@ngrad))) {
            cat("dimension of si is",dim(object$si),", but we want",c(object@ddim,object@ngrad),"\n")
            return(invisible(FALSE))
          }
          if (length(object@s0ind)<1) {
            cat("no S_0 images, parameters not identifiable \n")
             return(invisible(FALSE))
         }
         }
         )
setClass("dtiTensor",
         representation(method = "character"),
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
          if (!(object@method %in% c("linear","nonlinear"))) {
            cat("method should specify either linear or nonlinear \n")
            return(invisible(FALSE))
          }
          if (object@method=="nonlinear"&&(is.null(object$Varth)
               ||length(dim(object$Varth))!=4||
                 any(dim(object$Varth)[1:3]!=object@ddim)||
                 dim(object$Varth)[4]!=28)) {
            cat("no or incorrect component Varth in list \n")
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
