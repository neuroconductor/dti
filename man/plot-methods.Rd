\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY-method}
\alias{plot,dti-method}
\alias{plot,dtiData-method}
\alias{plot,dtiIndices-method}
\alias{plot,dtiTensor-method}
\alias{\method{plot},dtiData}
\alias{\method{plot},dtiIndices}
\alias{\method{plot},dtiTensor}
\title{2D Visualizing of dti-objects}
\usage{
\methods{plot}{dtiData}function(x, y,slice=1, gradient=NULL, view= "axial", show=TRUE,xind=NULL,yind=NULL,zind=NULL, mar=c(3,3,3,.3),mgp=c(2,1,0), ...)
\methods{plot}{dtiTensor}(x, y, slice=1, view="axial", quant=0, minanindex=NULL, contrast.enh=1,what="FA", qrange=c(.01,.99),xind=NULL,yind=NULL,zind=NULL, mar=c(2,2,2,.2),mgp=c(2,1,0),...)
\methods{plot}{dtiIndices}(x, y, slice=1, view= "axial", method=1, quant=0, minanindex=NULL, show=TRUE, contrast.enh=1,what="FA",xind=NULL,yind=NULL,zind=NULL, mar=c(3,3,3,.3),mgp=c(2,1,0), ...)
}
\arguments{
  \item{x}{ object of class \code{dtiIndices}, \code{dtiData} or \code{dtiTensor}}
  \item{y}{ not used }
   \item{slice }{Slice number}
   \item{view}{ choose sagittal, coronal, or axial view here}
   \item{gradient}{index of data cube to chose.} 
   \item{method }{\code{method==1} and \code{method==2} specify two different color schemes for 
    directional FA-maps. \code{method==3} specifies visualization of dtiIndices using color coded bary coordinates. }
   \item{quant }{if \code{is.null(minanindex)} specify \code{minanindex} as corresponding quantile of the fractal anisotropy (FA) index.}
   \item{minanindex }{display only information for voxel with \code{FA>minanindex}}
   \item{show }{Visualize information in a graphics device }
   \item{contrast.enh}{ enhance image contrast using \code{min(1,x$anindex/contrast.enh} instead
    of the anisotropy index itself. Effective values are within the interval (0,1).}
  \item{what}{ If \code{what="GA"} use geodesic anisotropy (GA) instead of 
   fractional anisotropy (FA).}
   \item{xind}{If provided restrict display to indices specified in \code{xind} for x-direction.}
   \item{yind}{If provided restrict display to indices specified in \code{xind} for y-direction.}
   \item{zind}{If provided restrict display to indices specified in \code{xind} for z-direction.}
   \item{\dots}{No further arguments implemented.}
}
\description{
 Methods for function \code{plot} in Package `graphics' 
}
\section{Methods}{
\describe{

\item{x = "ANY"}{ generic plot function }

This function allows to select a slide (2D) by arguments 
\code{slice} and \code{view}. 
\item{x = "dti"}{ see specific class. }

\item{x = "dtiData"}{\code{gradient} can be used to specify a specific data cube
associated with the index of a gradient direction. Images are produces that are
scaled by the maximal observed image value. This guarantees 
that subsequently produced images are on a comparable grey scale.
The resulting image is returned.}

\item{x = "dtiIndices"}{ Color coded anisotropy
maps are produced depending on the specification in \code{method}.
The resulting image is returned.}

\item{x = "dtiTensor"}{The tensor itself, fractional anisotropy (FA), mean diffusivity (MD) and a color coded anisotropy map are provided. \code{NULL} is returned.  }
}}
\value{
an image of class \code{adimpro} for objects of class \code{"dtiData"} or \code{"dtiIndices"}, \code{invisible(NULL)} otherwise.
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiIndices}}, \code{\link{dtiData}}, \code{\link{dtiTensor}}}
\examples{
\dontrun{demo(dti_art)}
}
\keyword{methods}
\keyword{hplot}
