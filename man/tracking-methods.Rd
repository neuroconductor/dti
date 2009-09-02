\name{tracking-methods}
\docType{methods}
\alias{tracking}
\alias{tracking-methods}
\alias{tracking,ANY-method}
\alias{tracking,dtiTensor-method}
\alias{tracking,dtiIndices-method}
\title{ Methods for Function `tracking' in Package `dti' }
\description{
  The function provides fiber tracking and 3D visualization of \code{"dtiTensor"} objects using the \code{"rgl"}-package. Functionality of the rgl-package allows to rotate and zoom the displayed object.
}
\usage{
  \S4method{tracking}{dtiTensor}(obj, xind=NULL, yind=NULL, zind=NULL, method="LINEPROP", minanindex=0.3, maxangle=30)
}
\arguments{
\item{obj}{ An object of class \code{dtiTensor}}
\item{xind}{Indices defining the ROI in x direction. Currently min/max is used to define ROIx}
\item{yind}{Indices defining the ROI in y direction. Currently min/max is used to define ROIy}
\item{zind}{Indices defining the ROI in z direction. Currently min/max is used to define ROIz}
\item{method}{Method for fibre tracking. "LINEPROP" is simple line propagation algorithm which is the default.}
\item{minanindex}{Minimal FA to follow the tracks. default 0.3}
\item{maxangle}{Maximal angle between fiber in adjacent voxels. default 30 degree.}
}
\section{Methods}{
\describe{

  \item{obj = "dtiTensor"}{ Fiber tracking is performed on the estimated vector field of principal diffusion direction using the method \code{method}. Currently only line propagation is implemented. The resulting tracks can be visualized using function \code{\link{show3d}}}. 
 
  \item{obj = "dtiIndices"}{ Fiber tracking is performed on the estimated vector field of principal diffusion direction using the method \code{method}. Currently only line propagation is implemented. The resulting tracks can be visualized using function \code{\link{show3d}}}. 
  Displays can be closed using function \code{rgl.close}
}
}
\value{
  The function returns an object of class \code{\linkS4class{dwiFiber}}.
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}, Joerg Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{dtiTensor}},\code{\linkS4class{dtiIndices}},\code{\linkS4class{dwiFiber}},
  \code{\link{show3d}}, \code{\link{summary}}, \code{\link{print}}
}
\keyword{ methods }
