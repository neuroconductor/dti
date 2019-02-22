\name{AdjacencyMatrix}
\alias{AdjacencyMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Create an adjacency matrix from fiber tracking results
}
\description{
   The function takes two objects, \code{fiberobj} with class 'dwiFiber' containing
   fiber tracking results and an array or nifti-object containing atlas information.
   For each combination of regions defined in the atlas the number of fibers connecting these
   regions is calculated, resulting in a matrix of fiber counts. As default this matrix
   is standardized and the diagonal elements are set to zero.
}
\usage{
AdjacencyMatrix(fiberobj, atlas, labels = NULL,
  method = c("standardize", "counts"), diagelements = FALSE,
  symmetric=TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fiberobj}{
     an object of class 'dwiFiber'
}
  \item{atlas}{
     an object of class 'array' or 'nifti' containing region indices as
    intensities. The atlas needs to be registered to DWI (subject) space,
    with array dimension corresponding to  \code{fiberobj@ddim}
}
  \item{labels}{
    optional labels for the regions. Will be used as dimnames of the resulting matrix.
}
  \item{method}{
   either \code{"standardize"} or \code{"counts"}, determines if fiber counts
   or a standardized (default) matix is returned.
}
  \item{diagelements}{
   logical, if \code{FALSE} the diagonal elements of the standardized matrix are set to zero (default).
}
\item{symmetric}{
  logical, with \code{ni} the number of fibers originating if \code{FALSE} standardized values
  \code{counts(i,j)/ni},
  if \code{TRUE} we get \code{counts(i,j)/sqrt(nj*nj)}.
}
  \item{verbose}{
   logical, if \code{TRUE} report pairwise fiber counts.
}
}
\value{
A matrix with dimensions equal to the number of regions defined in the atlas
and dimnames given by \code{labels} or by the region number. The matrix contains
fiber counts or values standardized with the number of fibers \code{ni, nj}
 originating/ending from the pair of regions. Depending on \code{symmetric} standardization
 is with \code{1/sqrt(ni*nj)} or with \code{1/ni}.
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\linkS4class{dwiFiber}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }% __ONLY ONE__ keyword per line
