\name{dtiIndices-methods}
\docType{methods}
\alias{dtiIndices-methods}
\alias{dtiIndices,ANY-method}
\alias{dtiIndices,dtiTensor-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Methods for Function `dtiIndices' in Package `dti' }
\description{
  The method creates estimates of the fractional anisotropy (FA) and relative anisotropy (RA) indices, the main directions of anisotropy and several statistics used for visualization.
}
\usage{
  \S4method{dtiIndices}{dtiTensor}(object, which)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object of class \code{dtiTensor} }
  \item{which}{ which indices should be created, currently not implemented. }
}
\value{an object of class \code{dtiIndices}.
}
\section{Methods}{
\describe{

\item{obj = "ANY"}{ Returns a warning. }

\item{obj = "dtiTensor"}{ Estimate tensor indices like trace, fractional and geodesic anisotropy, main diffusion direction and shape parameters. }

}
}
\references{  K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).

              \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiTensor-methods}}, \code{\linkS4class{dtiTensor}}, \code{\linkS4class{dtiIndices}}}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ manip }
\keyword{ models }
