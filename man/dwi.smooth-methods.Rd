\name{dwi.smooth-methods}
\docType{methods}
\alias{dwi.smooth}
\alias{dwi.smooth.ms}
\alias{dwi.smooth-methods}
\alias{dwi.smooth,ANY-method}
\alias{dwi.smooth.ms,ANY-method}
\alias{dwi.smooth,dtiData-method}
\alias{dwi.smooth.ms,dtiData-method}
\title{Smooth DWI data}
\description{
 Adaptive smoothing of DWI data. Smoothing is performed both in space and on the sphere (e.g. between images  obtained for different gradient directions) employing a natural
geometrical distance ( in SE(3)). Structural adaptation is used in space only.
Method \code{dwi.smooth} refers to the original POAS approach for single shell data.
Method \code{dwi.smooth.ms} implements an improved method that is applicable for both single and multi-shell data.
}
\usage{
  \S4method{dwi.smooth}{dtiData}(object, kstar, lambda=20, kappa0=NULL, mask=NULL, ncoils=1,
                                 sigma=NULL, level=NULL, vred=4, verbose=FALSE, dist=1, 
                                 model=c("Gapprox","Gapprox2","Chi","Chi2"))
  \S4method{dwi.smooth.ms}{dtiData}(object, kstar, lambda=12, kappa0=.5, ncoils=1,
                                    sigma=NULL, ws0=1, level=NULL, mask = NULL, xind=NULL,
                                    yind=NULL, zind=NULL, verbose=FALSE,
                                    usemaxni=TRUE, memrelease = TRUE)
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{kstar}{ Number of steps in structural adaptation}
  \item{lambda}{Scale parameter in adaptation }
  \item{kappa0}{determines amount of smoothing on the sphere. Larger values correspond to stronger smoothing on the sphere. If \code{kappa0=NULL} a value is that corresponds to a variace reduction with factor \code{vred} on the sphere.}
  \item{ncoils}{Number of coils in MR system }
  \item{sigma}{Error standard deviation. Assumed to be known and homogeneous in the current implementation. A reasonable estimate may be defined
  as the modal value of standard deviations obtained using method \code{getsdofsb}. }
  \item{level}{Threshold for image intensities when setting mask.}
  \item{mask}{Binary 3D image defining a mask}
  \item{vred}{Used if \code{kappa0=NULL} to specify the variance reduction on the sphere when suggesting a value of kappa0}
  \item{xind}{index for x-coordinate}
  \item{yind}{index for y-coordinate}
  \item{zind}{index for z-coordinate}
  \item{verbose}{If \code{verbose=TRUE} additional reports are given.}
  \item{dist}{Distance in SE3. Reasonable values are 1 (default, see Becker et.al. 2012), 2 ( a slight modification of 1: with k6^2 instead of abs(k6)) and 3 (using a 'naive' distance on the sphere)}
  \item{model}{Determines which quantities are smoothed. Possible values are
  \code{"Chi"} for observed values (assumed to be distributed as noncentral
  Chi with \code{2*ncoils} degrees of freedom), \code{"Chi2"} for squares of
  observed values (assumed to be distributed as noncentral
  Chi-squared with \code{2*ncoils} degrees of freedom). \code{"Gapprox"}
  and \code{"Gapprox2"} use a Gaussian approximation for the noncentral
  Chi distribution to smooth ovserved and squared values, respectively.
  }
  \item{ws0}{Factor to downweight information from S0 images, defaults
   to \code{1/numer of s0 images}.}
  \item{usemaxni}{If \code{"usemaxni==TRUE"} a strikter penalization is used.}
   \item{memrelease}{If \code{"memrelease==TRUE"} try to release allocated memory whenever possible.}
  }
\value{
  An object of class \code{"dtiData"} with smoothed diffusion weighted images.
}
\section{Methods}{
\describe{

\item{\code{signature(object) = "ANY"}}{
Returns a warning.
}

\item{\code{signature(object) = "dtiData"}}{
Smoothing of DWI data
}
}}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
  \code{\link{dtiData}},
  \code{\linkS4class{dtiData}},
}
\keyword{methods}
\keyword{smooth }
