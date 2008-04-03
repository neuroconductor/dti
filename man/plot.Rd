\name{plot}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualizing of dtiIndices}
\description{
Generic plot function for visualization of objects of class
\code{dtiIndices}. 
}
\usage{
plot(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{dtiIndices}, \code{dtiData} or \code{dtiTensor}}
  \item{y}{ not used }
  \item{\dots}{ additional arguments such as
   \item{slice }{Slice number}
   \item{view}{ choose sagittal, coronal, or axial view here}
   \item{method }{\code{method==1} and \code{method==2} specify two different color schemes for 
    directional FA-maps. \code{method==3} specifies visualization of dtiIndices using color coded bary coordinates}
   \item{quant }{if \code{is.null(minanindex)} specify \code{minanindex} as corresponding quantile of the fractal anisotropy (FA) index.}
   \item{minanindex }{display only information for voxel with \code{FA>minanindex}}
   \item{show }{Visualize information in a graphics device}
   \item{contrast.enh}{ enhance image contrast using \code{min(1,x$anindex/contrast.enh} instead
    of the anisotropy index itself. Effective values are within the interval (0,1).
 }
}
}
\details{
This function will be extended to visualize objects of class \code{dtiData}
and \code{dtiTensor}.
}
\value{
an image of class \code{adimpro}.
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\note{ }
\seealso{ \code{\link{dtiIndices}}, \code{\link{dtiData}}, \code{\link{dtiTensor}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
