\name{dwiRiceBias}
\docType{methods}
\alias{dwiRiceBias}
\alias{dwiRiceBias-methods}
\alias{dwiRiceBias,ANY-method}
\alias{dwiRiceBias,dtiData-method}
\title{Correction of Rician Bias for dwiData objects.
}
\description{
Correction of Rician Bias for dwiData objects. Noise standard deviation is assumed
to be known.
}
\usage{
\S4method{dwiRiceBias}{dtiData}(object, sigma=NULL, method="1stMoment")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{"dtiData"}}
  \item{sigma}{
  error standard deviation}
  \item{method}{
  method used for bias correction. Currently only 
  one method is implemented. It requires the error
  standard deviation to be known and the uses the 
  equation for the first moment of the Rician distribution.
}
}
\details{
  Bias correction based on known error standard deviation
  and first order moments.
}
\value{
   an object of class \code{"dtiData"} with bias corrected values
   in slot \code{si}
}
\section{Methods}{
\describe{

  \item{obj = "ANY"}{ Returns a warning. }

  \item{obj = "dtiData"}{ Returns an object of class \code{"dtiData"} with bias corrected values
   in slot \code{si}. }

}
}
\references{
             J. Polzehl and K. Tabelow, 
             \emph{Beyond the diffusion tensor model: The package \pkg{dti}},
             Journal of Statistical Software, to appear.

             K. Tabelow, H.U. Voss and J. Polzehl, 
             \emph{Modeling the orientation distribution function by 
             mixtures of angular central Gaussian distributions},
             Journal of Neuroscience Methods, to appear.

             J. Polzehl and K. Tabelow, 
             \emph{Structural adaptive smoothing in diffusion tensor imaging: 
             The R package dti}, Journal of Statistical Software, 31 (2009) pp. 1--24.
             
              K. Tabelow, J. Polzehl, V. Spokoiny and H.U. Voss. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).

             \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{dtiData}}, 
  \code{\link{readDWIdata}}, 
  \code{\link{dtiIndices-methods}},  
  \code{\link{medinria}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}}
  \code{\linkS4class{dwiMixTensor}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  \dontrun{demo(dti_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ manip }
\keyword{ models }
