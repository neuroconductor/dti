\name{clipdti-methods}
\docType{methods}
\alias{clipdti}
\alias{clipdti-methods}
\alias{clipdti,ANY-method}
\alias{clipdti,dtiData-method}
\alias{clipdti,dtiIndices-method}
\alias{clipdti,dtiTensor-method}
\title{Index operations on dti-objects}
\description{
 Generic methods for objects of class \code{dti}. It allows to
  extract part of the data contained in the object by specifying
  a subcube. A new object of the same class containing the reduced data
  is created.
}
\usage{
clipdti(obj, ...)
\method{clipdti}{dtiData}(obj,xind=NULL,yind=NULL,zind=NULL)
\method{clipdti}{dtiIndices}(obj,xind=NULL,yind=NULL,zind=NULL)
\method{clipdti}{dtiTensor}(obj,xind=NULL,yind=NULL,zind=NULL)
}
\arguments{
  \item{obj}{ an object of class \code{"dtiData"}, \code{"dtiTensor"} or
     \code{"dtiIndices"}}
  \item{xind}{ vector of x-coordinates, defaults to whole range. }
  \item{yind}{ vector of y-coordinates, defaults to whole range. }
  \item{zind}{ vector of z-coordinates, defaults to whole range. }
  \item{\dots}{ Currently no additional arguments are implemented.}
}
\section{Methods}{
\describe{

\item{obj = "ANY"}{ returns a warning. }

\item{obj = "dtiData"}{ Returns an valid object of class \code{dtiData} containing reduced data.}

\item{obj = "dtiIndices"}{ Returns an valid object of class \code{dtiIndices} containing reduced data. }

\item{obj = "dtiTensor"}{ Returns an valid object of class \code{dtiTensor} containing reduced data. }
Information is extracted for voxel 
   within the  cube defined by
   arguments \code{xind}, \code{yind} and \code{zind}.
}}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{\code{\link{dtiTensor-class}}, \code{\link{dtiIndices-class}}  }
\examples{
\dontrun{demo(dti_art)}
}
\keyword{methods}
\keyword{manip}
