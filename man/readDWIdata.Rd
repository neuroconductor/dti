\name{readDWIdata}
\alias{readDWIdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read Diffusion Weighted Data }
\description{
  The function reads Diffusion Weighted Data from medicial imaging file in data directories for further processing.
}
\usage{
readDWIdata(dirlist, format, nslice, gradient, order = NULL, xind = NULL, yind = NULL, zind = NULL, level = 0, mins0value = 0, maxvalue = 10000, voxelext = NULL, orientation = c(1, 3, 5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dirlist}{ list of directories containing the data files }
  \item{format}{ string specifying the medical imaging format, one of ''DICOM'', ''NIFTI'', ''ANALYZE'', or ''AFNI'' }
  \item{nslice}{ number of slices (usually z-direction) in dataset }
  \item{gradient}{ matrix of gradient vectors }
  \item{order}{ vector, specifying a different order of the data files, i.e. other than alphabetic order in the directories given by \code{dirlist}. If not given, 1:n is used for n data files (no order change). }
  \item{xind}{ subindex for x-direction }
  \item{yind}{ subindex for y-direction }
  \item{zind}{ subindex for z-direction }
  \item{level}{ determine \code{mins0value} as quantile of positive S0-values }
  \item{mins0value}{ set voxel in S0-images with values less than \code{level} inactive }
  \item{maxvalue}{ set voxel with values larger than \code{maxvalue} inactive }
  \item{voxelext}{ voxelextensions in coordinate directions }
  \item{orientation}{ orientations of data coded as in AFNI }
}
\details{
  The function reads the data files given in the directories in \code{dirlist} in alphabetic order. The order can be changed using the \code{order} argument: If \code{filelist} is the vector of files in alphabetic order, they are read in the order \code{filelist[order]}. If order is not given \code{order <- 1:n} is used (no change!). The medical imaging format is given by \code{format} and can be one of ''DICOM'', ''NIFTI'', ''ANALYZE'', or ''AFNI''. The number of slices of the three dimensional data cube is given by \code{nslice}. The diffusion gradients are provided as matrix \code{gradient}.

  \code{xind}, \code{yind}, and \code{zind} define a region of interest as indices. If not given \code{1:dim[i]} is used. \code{level} determine \code{mins0value} as quantile of positive S0-values. \code{mins0value} sets voxel in S0-images with values less than \code{level} inactive. \code{maxvalue} sets voxel with values larger than \code{maxvalue} inactive.

  \code{voxelext} defines the voxel extension, overwrites the values found in the imaging files. \code{orientation} codes the data orientation in AFNI notation.
}
\value{
  An object of class \code{dtiData}.
}
\references{  K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, NeuroImage 39(4), 1763-1773 (2008).

\url{http://www.wias-berlin.de/projects/matheon_a3/} }
\author{Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de}
\seealso{  \code{\link{dtiTensor}}, \code{\linkS4class{dtiData}} }
\examples{
\dontrun{demo(dti_art)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
