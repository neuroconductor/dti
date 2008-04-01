\name{dtiData}
\alias{dtiData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dtiData(gradient, imagefile, ddim, xind = NULL, yind = NULL, zind = NULL, level = 0, mins0value = 0, maxvalue = 10000, voxelext = c(1, 1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gradient}{ ~~Describe \code{gradient} here~~ }
  \item{imagefile}{ ~~Describe \code{imagefile} here~~ }
  \item{ddim}{ ~~Describe \code{ddim} here~~ }
  \item{xind}{ ~~Describe \code{xind} here~~ }
  \item{yind}{ ~~Describe \code{yind} here~~ }
  \item{zind}{ ~~Describe \code{zind} here~~ }
  \item{level}{ ~~Describe \code{level} here~~ }
  \item{mins0value}{ ~~Describe \code{mins0value} here~~ }
  \item{maxvalue}{ ~~Describe \code{maxvalue} here~~ }
  \item{voxelext}{ ~~Describe \code{voxelext} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/}}
\author{Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\note{}
\seealso{  \code{\link{dtiTensor}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
