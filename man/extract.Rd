\name{extract}
\alias{extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract information from objects of class dtiData, dtiTensor or dtiIndices }
\description{
  The function extracts or/and computes specified statistics from 
object of class \code{dti}. This can be restricted to a subset of voxel.
}
\usage{
extract(obj, ...)
\method{extract}{dtiData}(obj,what="dtiData",xind=NULL,yind=NULL,zind=NULL)
\method{extract}{dtiTensor}(obj,what="dtiData",xind=NULL,yind=NULL,zind=NULL)
\method{extract}{dtiIndices}(obj,what="dtiData",xind=NULL,yind=NULL,zind=NULL)
}
%- maybe also 'usage' for other objects documented here.
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
\details{
   Extraction of statistics is restricted to the cube defined by
   arguments \code{xind}, \code{yind} and \code{zind}.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
