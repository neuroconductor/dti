\name{clipdti}
\alias{clipdti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Index operations on dti-objects }
\description{
  Generic function for objects of class \code{dti}. It allows to
  extract part of the data contained in the object by specifying
  a subcube. A new object of the same class containing the reduced data
  is created.
}
\usage{
clipdti(obj, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ an object of class \code{"dtiData"}, \code{"dtiTensor"} or
     \code{"dtiIndices"}}
  \item{\dots}{ Currently implemented are 
  arguments \code{xind}, \code{yind}, \code{zind}. All arguments default to 
  \code{NULL} and allow to specify an index vector for the respective dimension}
}
\details{
3D data volumes are reduced by applying the index vectors 
supplied in \code{xind}, \code{yind}, \code{zind}. If an index vector
is NULL then now reduction is performed for the corresponding dimension.
An new object is created that inherits the class from \code{obj}.
}
\value{
An object of class \code{class(obj)}. Slots \code{xind}, \code{yind} and \code{zind} contain the indices with respect to the original object.
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).
 }
\url{http://www.wias-berlin.de/projects/matheon_a3/} 
\author{Karsten Tabelow tabelow@wias-berlin.de, Joerg Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiData}}, \code{\link{dtiTensor}}, \code{\link{dtiIndices}}, \code{\link{dtiData-class}}, \code{\link{dtiTensor-class}}, \code{\link{dtiIndices-class}}}
\examples{
\dontrun{demo(dti_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
