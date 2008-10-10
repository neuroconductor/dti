\name{print-methods}
\docType{methods}
\alias{print}
\alias{print-methods}
\alias{print,ANY-method}
\alias{print,dti-method}
\title{ Print Method }
\description{
 The function provides information on data dimensions, data source and
 existing slot-names for objects of class \code{dti}, \code{dtiData},
 \code{dtiTensor} and \code{dtiIndices}.
}
\usage{
print(x, ...)
\S4method{print}{dti}(x)
}
\arguments{
  \item{x}{ object of class \code{dtiIndices}, \code{dtiData} or \code{dtiTensor}}
  \item{\dots}{ additional arguments}
}
\section{Methods}{
\describe{

\item{x = "ANY"}{ Returns a warning }

\item{x = "dti"}{ The function provides information on data dimensions,
  data source and existing slot-names for objects of class \code{dti}. }
}}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\linkS4class{dtiIndices}}, \code{\linkS4class{dtiData}}, \code{\linkS4class{dtiTensor}}}
\examples{
\dontrun{demo(dti_art)}
}
\keyword{methods}
\keyword{ utilities}
