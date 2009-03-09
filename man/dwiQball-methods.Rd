\name{dwiQball-methods}
\docType{methods}
\alias{dwiQball}
\alias{dwiQball-methods}
\alias{dwiQball,ANY-method}
\alias{dwiQball,dtiData-method}
\title{ Methods for Function `dwiQball' in Package `dti' }
\description{ 
  The method estimates, in each voxel, the coefficients of an expansion of the apparent diffusion cefficient (ADC) with respect to a apherical harmonics  orthonormal system  from the DWI data contained in an object of class \code{"dtiData"}.
}
\usage{
  \S4method{dwiQball}{dtiData}(object, method="linear", order=0, lambda=0)
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{method}{ Method for tensor estimation. May be \code{"linear"}, or \code{"nonlinear"}. }
  \item{order}{ even integer: maximum order of the sperical harmonics expansion }
  \item{lambda}{ nonnegative regularization parameter. }
}
\value{ 
  An object of class \code{"dwiQball"}.
}
\section{Methods}{
\describe{

  \item{obj = "ANY"}{ Returns a warning. }

  \item{obj = "dtiData"}{ Estimate, in each voxel, the coefficients of an expansion of the apparent diffusion cefficient (ADC) with respect to a apherical harmonics  orthonormal system.  If \code{method=="linear"} estimates are obtained using a linearization of the ADC model, see Descoteux (2007). \code{method=="nonlinear"}
not yet implemented.  }

}
}
\references{ M. Descoteaux, E. Angelino, S. Fitzgibbons and R. Deriche,
  \emph{Regularized, Fast and Robust Analytical Q-Ball Imaging},
  Magnetic Resonance Methods, 2007, 58, 497-512.
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
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
\keyword{ manip }
\keyword{ models }
