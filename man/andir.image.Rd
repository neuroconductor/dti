\name{andir.image}
\alias{andir.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Illustrate anisotropy within DT images.}
\description{
  The function generates an image characterizing the anisotropy within a specified
slide of an 3D diffusion tensor image. 
}
\usage{
andir.image(dtobject, slice = 1, method = 1, quant = 0, minanindex = NULL, show=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtobject}{ object of class \code{dti} as produced by function \code{\link{dtianiso}}. }
  \item{slice}{ slice number }
  \item{method}{ \code{method} specifies the color coding of the anisotropy direction. }
  \item{quant}{ If \code{is.null(minanindex)} then \code{quant} 
       determines the lower bound for the anisotropy index as the \code{quant}-quantile of estimated
       anisotropy indices.  }
  \item{minanindex}{\code{minanindex} allows to determine the lower bound for the anisotropy index. 
   only results at voxel with anisotropy index larger than this lower bound will be displayed. }
  \item{show}{ Determines if the resulting object is displayed.}
  \item{...}{ Additional parameters for function \code{show.image} from package \code{adimpro} or function \code{image}.}
}
\details{
}
\value{
  If package \code{adimpro} is installed the function creates an image object of class \code{"adimpro"}, otherwise
  an array of rgb-values, characterizing the anisotropy within the specified slice. 
}
\references{}
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de }
\note{ 
}
\seealso{\code{\link{dtianiso}} }
\examples{
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ manip }% __ONLY ONE__ keyword per line
