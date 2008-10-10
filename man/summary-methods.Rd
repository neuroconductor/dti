\name{summary-methods}
\docType{methods}
\alias{summary}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{summary,dti-method}
\title{ Method for summary }
\description{
  The method provides summary information for objects of class \code{dti}.
}
\usage{
summary(object, ...)
\S4method{summary}{dti}(object)
}
\arguments{
  \item{object}{ an object of class \code{dti}, \code{dtiData},
\code{dtiTensor} or \code{dtiIndices}. }
  \item{\dots}{ without effect }
}
\section{Methods}{
\describe{

\item{object = "ANY"}{ Returns a warning }

\item{object = "dti"}{  The function provides summary information for objects of class \code{dti}, \code{dtiData},
\code{dtiTensor} and \code{dtiIndices}. }
}}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiIndices}}, \code{\link{dtiData}}, \code{\link{dtiTensor}}}
\examples{
\dontrun{demo(dti_art)}
}
\keyword{methods}
\keyword{utilities}
