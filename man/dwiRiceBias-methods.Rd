\name{dwiRiceBias-methods}
\docType{methods}
\alias{dwiRiceBias-methods}
\alias{dwiRiceBias,ANY-method}
\alias{dwiRiceBias,dtiData-method}
\title{ Correction for Rician Bias }
\description{
 Correction for Rician Bias assuming known variance parameter
}
\usage{
  \S4method{dwiRiceBias}{dtiData}(object, sigma=NULL, method="1stMoment")
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{sigma}{ Value for the second parameter of the Rician distribution. }
  \item{method}{ method used to estimate the first parameter of the Rician distribution.}
}
\value{ 
  An object of class \code{"dtiData"}.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{ Returns a warning. 
}

\item{\code{signature(object = "dtiData")}}{
Returns a \code{dtiData} object with bias-corrected image intensities.
}
}}
\keyword{methods}
\keyword{ manip }
