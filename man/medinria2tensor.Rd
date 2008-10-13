\name{medinria2tensor}
\alias{medinria2tensor}
\title{ Read Diffusion Tensor Data from NIFTI File }
\description{
  Read Tensor-Data from NIFTI File. Interface function to MedINRIA 
}
\usage{
medinria2tensor(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ name of the file containing the tensor data. }
}
\details{
The function allows to read tensor data written as NIFTY-Files by MedINRIA.
}
\value{
  Object of class \code{dtiTensor}.
}
\references{  K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).

              \url{ http://www.wias-berlin.de/projects/matheon_a3/}
 
              P. Fillard, J. Souplet and N. Toussaint 
              \emph{ Medical Image Navigation and Research Tool by INRIA (MedINRIA)}, 
              INRIA Sophia Antipolis - Research Project ASCLEPIOS 2007
}
\author{Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de}
\seealso{ \code{\link{dtiTensor-methods}} }
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ IO }
