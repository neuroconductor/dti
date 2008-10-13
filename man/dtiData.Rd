\name{dtiData}
\alias{dtiData}
\title{ Read Diffusion Weighted Data from Image File}
\description{
  This function creates a \code{"dtiData"} object from an imagefile, where the diffusion weighted data is given as 2-byte integer.
}
\usage{
  dtiData(gradient, imagefile, ddim, xind = NULL, yind = NULL, zind = NULL, level = 0, mins0value = 0, maxvalue = 10000, voxelext = c(1, 1, 1), orientation = c(1, 3, 5))
}
\arguments{
  \item{gradient}{ matrix of diffusion gradients (including zero gradients for S0 images}
  \item{imagefile}{ filename for image data (binary 2Byte integers)}
  \item{ddim}{ dimension of image cube (3D) }
  \item{xind}{ subindex for x-direction }
  \item{yind}{ subindex for y-direction }
  \item{zind}{ subindex for z-direction }
  \item{level}{ determine \code{mins0value} as quantile of positive S0-values  }
  \item{mins0value}{ set voxel in S0-images with values less than \code{level} inactive }
  \item{maxvalue}{ set voxel with values larger than \code{maxvalue} inactive  }
  \item{voxelext}{ voxelextensions in coordinate directions.}
  \item{orientation}{ orientations of data coded as in AFNI}
}
\details{
  This function creates a dtiData object from an imagefile, where the diffusion weighted data is given as 2-byte integer. This image file has to be prepared by the user and there is not yet a function to do this within the package. Simply use \code{\link{writeBin}} to write out first all S0 images and than all Si images. The \code{gradient} should be created according to this order. Run the demo in order to have an example, how to do this!

  \code{xind}, \code{yind}, and \code{zind} define a region of interest as indices. If not given \code{1:dim[i]} is used. \code{level} determine \code{mins0value} as quantile of positive S0-values. \code{mins0value} sets voxel in S0-images with values less than \code{level} inactive. \code{maxvalue} sets voxel with values larger than \code{maxvalue} inactive.

  \code{voxelext} defines the voxel extension. 
}
\value{
An object of class \code{dtiData}.
}
\references{  K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).

              \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{ \code{\link{dtiTensor-methods}}, \code{\linkS4class{dtiData}} }
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ IO }
