\name{show3d-methods}
\docType{methods}
\alias{show3d-methods}
\alias{show3d,ANY-method}
\alias{show3d,dtiIndices-method}
\alias{show3d,dtiTensor-method}
\title{ Methods for Function `show3d' in Package `dti' }
\description{
  The function provides 3D visualization of \code{"dtiTensor"} and \code{"dtiIndices"} objects using the \code{"rgl"}-package. Functionality of the rgl-package allows to rotate and zoom the displayed object.
}
\usage{
  \S4method{show3d}{dtiTensor}(obj, nx=NULL, ny=NULL, nz=NULL, center=NULL, method=1, level=0, scale=.25, bgcolor="black", add=FALSE, subdivide=2, smooth=TRUE, maxobjects=729, \dots)
  \S4method{show3d}{dtiIndices}(obj, index="FA", nx=NULL, ny=NULL, nz=NULL, center=NULL, method=1, level=0, bgcolor="black", add=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ An object of class \code{dtiTensor} or \code{dtiIndices} }
\item{nx}{Number of voxel in x-direction}
\item{ny}{Number of voxel in y-direction}
\item{nz}{Number of slices}
\item{center}{Vector of length 3 specifying the center of the data cube (class \code{dtiTensor}) or center of display (class \code{dtiIndices})}
\item{method}{\code{method==1} and \code{method==2} specify two different color schemes for 
    directional FA-maps.}
\item{level}{Minimal FA value of tensors.}
\item{scale}{Scale factor for the size of ellipsoids}
\item{bgcolor}{Backgroundcolor for rgl-display}
\item{add}{If true information is added to the current device, otherwise a new device is opened.}
\item{subdivide}{Level of subdivisions for meshing of ellipsoids, see package rgl.}
\item{smooth}{Use smooted surfaces, see package rgl.}
\item{maxobjects}{Maximal size of data cube (in voxel) to display}
\item{index}{Eiter \code{"FA"} for fractional anisotropy index or 
\code{"GA"} for geodetic anisotropy index.}
\item{\dots}{Additional parameters passed to functions from the rgl-package.}
}
\section{Methods}{
\describe{

\item{obj = "ANY"}{ Returns a warning }

\item{obj = "dtiIndices"}{ Objects of class \code{dtiIndices} are visualized as a collection of line segments with location given by the voxel center, orientation and color determined by the main direction of inisotropy and length corresponding to either fractional or geodesic anisotropy as specified in \code{index}. 

Displayed objects are restricted to voxel with an fractional (geodesic) anisotropy larger than \code{level}.
}

\item{obj = "dtiTensor"}{For objects of class \code{dtiTensor} ellipsoids are visualized at the voxel centers. Orientation and size correspond to the tensor
  values, color is determined by the main direction of anisotropy using the colorsceme specified with \code{method}. The fractional anisotropy value is coded as transparency. The value of \code{maxobjects} limits the size of datacube and may be increased
 on hardware with suitable graphics capabilities.
 }
Displays can be closed using function \code{rgl.close}
}
}
\value{
  The function returns the number of the current rgl-device.
}
\references{  K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).

              \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\linkS4class{dtiTensor}}, \code{\linkS4class{dtiIndices}}  }
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
\keyword{ iplot }
\keyword{ hplot }
