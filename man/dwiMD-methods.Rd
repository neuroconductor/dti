\name{dwiMD-methods}
\docType{methods}
\alias{dwiMD}
\alias{dwiMD-methods}
\alias{dwiMD,ANY-method}
\alias{dwiMD,dtiData-method}
\alias{dwiMD,dtiTensor-method}
\title{ Methods for Function `dwiMD' in Package `dti' }
\description{ 
  The method calculates, in each voxel, the mean diffusivityfrom the DWI data contained in an object of class \code{"dtiData"} or from the estimated diffusion tensor in an object of class \code{"dtiThensor"}
}
\usage{
  \S4method{dwiMD}{dtiData}(object, eps=.05)
  \S4method{dwiMD}{dtiTensor}(object)
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"} or class \code{"dtiTensor"} }
  \item{eps}{ In case of an \code{"dtiData"} object the method searches for three approximately
            orthogonal gradient directions. The mean diffusivity is then calculated using 
            the observed apparent diffusion coefficients for these gradients.
            The search is terminated if a measure for deviation from orthogonality is smaller than \code{eps}}
}
\value{ 
  An array of mean diffusivity values.
}
\section{Methods}{
\describe{

  \item{obj = "ANY"}{ Returns a warning. }

  \item{obj = "dtiData"}{ calculates mean diffusivity from 3 approximately orthogonal 
                          gradient directions}

  \item{obj = "dtiTensor"}{ calculates mean diffusivity from estimated diffusion tensors}
}
}
\references{ J. Polzehl and K. Tabelow, 
             \emph{Beyond the diffusion tensor model: The package \pkg{dti}},
             Journal of Statistical Software, 44(12), 1-26 (2011).

             J. Polzehl and K. Tabelow, 
             \emph{Structural adaptive smoothing in diffusion tensor imaging: 
             The R package dti}, Journal of Statistical Software, 31(9) 1-24 (2009).
             

              \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{dtiData}}, 
  \code{\link{readDWIdata}}, 
  \code{\link{dtiIndices-methods}},  
  \code{\link{medinria}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}}
  \code{\linkS4class{dwiMixtensor}}
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
\keyword{ manip }
\keyword{ models }
