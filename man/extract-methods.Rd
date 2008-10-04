\name{extract-methods}
\docType{methods}
\alias{extract}
\alias{extract-methods}
\alias{extract,ANY-method}
\alias{extract,dtiData-method}
\alias{extract,dtiIndices-method}
\alias{extract,dtiTensor-method}
\title{ Extract information from objects}
\description{
 The methods extract or/and compute specified statistics from 
object of class \code{dti}. This can be restricted to a subset of voxel.
}
\usage{
extract(obj, ...)
\method{extract}{dtiData}(obj,what="data",xind=NULL,yind=NULL,zind=NULL)
\method{extract}{dtiTensor}(obj,what="tensor",xind=NULL,yind=NULL,zind=NULL)
\method{extract}{dtiIndices}(obj,what=c("fa","andir"),xind=NULL,yind=NULL,zind=NULL)
}
\arguments{
  \item{obj}{ Object of class \code{dti} }
  \item{what}{ Statistic to extract. Options are

for \code{class(obj)=="dtiData"}
     \item{"btb"}{Gradient information}
     \item{"S0"}{S0 images}  
     \item{"Si"}{Sb images}  
     \item{"Sb"}{Sb images}  
     \item{"data"}{all images}

     for \code{class(obj)=="dtiTensor"}
     \item{"fa"}{Fractional anisotropy}  
     \item{"ga"}{Geodetic anisotropy}  
     \item{"md"}{Mean diffusivity}  
     \item{"evalues"}{Eigenvalues}
     \item{"andir"}{Main direction of anisotropy}
     \item{"tensor"}{Estimated tensor}
     \item{"s0"}{Estimated S0 image}
     \item{"mask"}{mask}

     for \code{class(obj)=="dtiIndices"}
     
     \item{"fa"}{Fractional anisotropy}  
     \item{"ga"}{Geodetic anisotropy}  
     \item{"md"}{Mean diffusivity}  
     \item{"andir"}{Main direction of anisotropy}
     \item{"bary"}{Bary coordinates}

A vector of options can be supplied.
}
  \item{xind}{ vector of x-coordinates, defaults to whole range. }
  \item{yind}{ vector of y-coordinates, defaults to whole range. }
  \item{zind}{ vector of z-coordinates, defaults to whole range. }
  \item{\dots}{ one or several  of the parameters above }
}
\section{Methods}{
\describe{

\item{obj = "ANY"}{ Returns a warning. }

\item{obj = "dtiData"}{    Extraction of squared gradient matrix or of
   S0, Sb or all images restricted to the cube defined by
   arguments \code{xind}, \code{yind} and \code{zind}.
 }

\item{obj = "dtiIndices"}{ Returns an array containing the specified statistics,
  i.e. fractional anisotropy (fa), geodesic anisotropy (ga), mean diffusivity (md),
    main direction of anisotropy (andir) and/or Bary coordinates (bary),
   as specified in argument \code{what}. Information is extracted for voxel 
   within the  cube defined by
   arguments \code{xind}, \code{yind} and \code{zind}.
 }

\item{obj = "dtiTensor"}{    Returns an array containing the specified statistics,
  i.e. fractional anisotropy (fa), geodesic anisotropy (ga), mean diffusivity (md),
   eigenvalues (evalues), main direction of anisotropy (andir), the tensor 
   (tensor) the estimated S0 image (s0) and/or the mask used to restrict computations (mask),
   as specified in argument \code{what}. Information is extracted for voxel 
   within the  cube defined by
   arguments \code{xind}, \code{yind} and \code{zind}.
 }
}
}
\value{
  List with components carrying the names of the options specified in
  argument \code{what}.
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} }
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{\code{\link{dtiTensor-class}}, \code{\link{dtiIndices-class}}  }
\examples{
\dontrun{demo(dti_art)}
}
\keyword{methods}
\keyword{manip}
