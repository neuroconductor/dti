\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY-method}
\alias{plot,dti-method}
\alias{plot,dtiData-method}
\alias{plot,dtiIndices-method}
\alias{plot,dtiTensor-method}
\title{ Methods for Function `plot' in Package `dti' } 
\description{
  Visualization of objects of class \code{"dtiData"}, \code{"dtiIndices"} and class \code{"dtiTensor"}
}
\usage{
  \S4method{plot}{dtiData}(x, y, slice=1, gradient=NULL, view= "axial", show=TRUE, xind=NULL, yind=NULL, zind=NULL, mar=c(3,3,3,.3), mgp=c(2,1,0), \dots)
  \S4method{plot}{dtiTensor}(x, y, slice=1, view= "axial", method=1, quant=0, minanindex=NULL, show=TRUE, contrast.enh=1, what="FA", xind=NULL, yind=NULL, zind=NULL, mar=c(3,3,3,.3), mgp=c(2,1,0), \dots)
  \S4method{plot}{dtiIndices}(x, y, slice=1, view="axial", quant=0, minanindex=NULL, contrast.enh=1,what="FA", qrange=c(.01,.99), xind=NULL, yind=NULL, zind=NULL, mar=c(2,2,2,.2), mgp=c(2,1,0), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{dtiIndices}, \code{dtiData} or \code{dtiTensor}}
  \item{y}{ not used }
  \item{\dots}{ additional arguments such as }
  \item{slice}{ Slice number }
  \item{view}{ choose sagittal, coronal, or axial view here }
  \item{gradient}{ index of data cube to chose. (for class \code{dtiData}) only. Defaults to the first S0 image.}
  \item{method}{ \code{method==1} and \code{method==2} specify two different color schemes for  directional FA-maps. \code{method==3} specifies visualization of dtiIndices using color coded bary coordinates. (for class \code{dtiIndices} only)}
  \item{quant}{ if \code{is.null(minanindex)} specify \code{minanindex} as corresponding quantile of the fractal anisotropy (FA) index.}
  \item{minanindex}{ display only information for voxel with \code{FA>minanindex}}
  \item{show}{Visualize information in a graphics device (for classes \code{dtiData} and \code{dtiIndices} only)}
  \item{contrast.enh}{ enhance image contrast using \code{min(1,x$anindex/contrast.enh} instead of the anisotropy index itself. Effective values are within the interval (0,1).}
  \item{what}{ If \code{what="GA"} use geodesic anisotropy (GA) instead of fractional anisotropy (FA).}
  \item{mar}{}
  \item{mgp}{}
  \item{qrange}{}
  \item{xind}{If provided restrict display to indices specified in \code{xind} for x-direction.}
  \item{yind}{If provided restrict display to indices specified in \code{xind} for y-direction.}
  \item{zind}{If provided restrict display to indices specified in \code{xind} for z-direction.} 
}
\value{ an object of class \code{dtiTensor}.
}
\section{Methods}{
\describe{
  \item{x = "ANY"}{ Generic function: see \code{\link[graphics]{plot}}. }

  \item{x =  "dti"}{ Returns a warning. }

  \item{x = "dtiData"}{ For objects of class \code{"dtiData"}
    \code{gradient} can be used to specify a specific data cube
    associated with the index of a gradient direction.
    For objects of class \code{"dtiData"} images are produces that are
    scaled by the maximal observed image value. This guarantees 
    that subsequently produced images are on a comparable grey scale.
    The resulting image is returned. }

  \item{x = "dtiIndices"}{ For objects of class \code{"dtiIndices"} color
    coded anisotropy maps are produced depending on the specification in
    \code{method}. The resulting image is returned. }

  \item{x = "dtiTensor"}{ For objects of class \code{"dtiTensor"} the
    tensor itself, fractional anisotropy (FA), mean diffusivity (MD) and a
    color coded anisotropy map are provided. \code{NULL} is returned. }
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
\keyword{ hplot }
