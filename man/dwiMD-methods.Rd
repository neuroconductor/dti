\name{dwiMD-methods}
\docType{methods}
\alias{dwiMD}
\alias{dwiMD,ANY-method}
\alias{dwiMD,dtiData-method}
\alias{dwiMD,dtiTensor-method}
\title{ Methods for Mean Diffusivity in Package `dti' }
\description{ Compute mean diffusivity (MD) from dtiData or dtiTensor objects
}
\usage{
  \S4method{dwiMD}{dtiData}(object, eps=.05)
  \S4method{dwiMD}{dtiTensor}(object)
}
\arguments{
  \item{object}{Object of class \code{"dtiData"} or \code{"dtiTensor"}}
  \item{eps}{tolerance in search for good gradient combinations.}
}
\value{
  Array of mean diffusivities.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{Returns a warning}

\item{\code{signature(object = "dtiData")}}{
   searches for three gradients that enable best MD evaluation. Returns MD.
}
\item{\code{signature(object = "dtiTensor")}}{
   calculates MD values from estimated difusion tensors
}
}
}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}
}
\seealso{
  \code{\link{dtiData}},
  \code{\link{dtiTensor}},
}
\keyword{methods}
\keyword{model}
