\name{dtiIndices-class}
\docType{class}
\alias{dtiIndices-class}
\alias{dtiTensor,dtiIndices-method}
\alias{plot,dtiIndices-method}

\title{Class "dtiIndices" }
\description{ Fractional Anisotropy indices and main directions of anisotropy derived from a diffusion tensor	}
\section{Objects from the Class}{
Objects can be created by calls of function \code{\link{dtiIndices}}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} usually empty }
    \item{\code{fa}:}{Object of class \code{"array"} Fractional anisotropy values (FA)}
    \item{\code{ra}:}{Object of class \code{"array"} Relative anisotropy  values (RA) }
    \item{\code{trc}:}{Object of class \code{"array"} mean eigenvalues of the diffusion tensor.  }
    \item{\code{bary}:}{Object of class \code{"array"} barycentric coordinates }
    \item{\code{lambda}:}{Object of class \code{"array"} eigenvalues of the diffusion tensor.}
    \item{\code{eigenv}:}{Object of class \code{"array"} eigenvectors of the diffusion tensor. }
\item{\code{btb}:}{Object of class \code{"matrix"} matrix of dimension \code{c(6,ngrad)} obtained from gradient directions.}
    \item{\code{ngrad}:}{Object of class \code{"integer"} number of gradients (including zero gradients) }
    \item{\code{s0ind}:}{Object of class \code{"integer"} index of zero gradients within sequence \code{1:ngrad} }
    \item{\code{replind}:}{Object of class \code{"integer"} index (identifier) of unique 
           gradient directions. Used to charactreize replications in the gradient design by identical indices. length is \code{ngrad}}
    \item{\code{ddim}:}{Object of class \code{"integer"} dimension of original image cubes. Integer vector of length 3. }
    \item{\code{ddim0}:}{Object of class \code{"integer"} dimension of subcube defined by \code{xind}, \code{yind} and \code{zind}. }
    \item{\code{xind}:}{Object of class \code{"integer"} index for subcube definition in x-direction  }
    \item{\code{yind}:}{Object of class \code{"integer"} index for subcube definition in y-direction }
    \item{\code{zind}:}{Object of class \code{"integer"} index for subcube definition in z-direction }
    \item{\code{voxelext}:}{Object of class \code{"numeric"} voxel extensions in x-, y- and z-direction. vector of length 3. }
    \item{\code{source}:}{Object of class \code{"character"} name of the imgfile used to create the data. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{dti}"}, directly.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "dtiIndices")}: 2D Visualization for a specified slice.  }
	 }
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} 
 }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de  }
\seealso{
\code{\link{dtiData}},\code{\link{dtiTensor}},\code{\link{dti.smooth}}, \code{\linkS4class{dti}}, \code{\linkS4class{dtiData}}, \code{\linkS4class{dtiTensor}}	
}
\examples{
showClass("dtiIndices")
}
\keyword{classes}
