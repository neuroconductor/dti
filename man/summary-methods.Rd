\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{summary,dti-method}
\title{ Methods for Function `summary' in Package `dti' }
\description{
  The method provides summary information for objects of class \code{dti}.
}
\usage{
  \S4method{summary}{dti}(object, ...)
}
\arguments{
  \item{object}{ an object of class \code{dti}, \code{dtiData},
\code{dtiTensor} or \code{dtiIndices}. }
  \item{\dots}{ Additional arguments in \code{\dots} 
are passed to function \code{quantile}, e.g. argument \code{probs} may be specified here.}
}
\section{Methods}{
\describe{

\item{object = "ANY"}{ Generic function: see \code{\link[base]{summary}}. }

\item{object = "dti"}{  The function provides summary information for objects of class \code{dti}, \code{dtiData},
\code{dtiTensor} and \code{dtiIndices}. }
}}
\references{  K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).

              \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\linkS4class{dtiIndices}}, \code{\linkS4class{dtiData}}, \code{\linkS4class{dtiTensor}}}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
\keyword{ utilities }
