\name{dti-package}
\alias{dti-package}
\alias{dti}
\docType{package}
\title{ DTI Analysis }
\description{
Diffusion Weighted Imaging is a Magnetic Resonance Imaging
modality, that measures diffusion of water in tissues like the human brain.
The package contains R-functions to process diffusion-weighted data in the
context of the diffusion tensor model (DTI). This includes the calculation
of anisotropy measures and, most important, the implementation of our structural
adaptive smoothing algorithm as described in K. Tabelow, J. Polzehl, V. Spokoiny,
and H.U. Voss, Diffusion Tensor Imaging: Structural Adaptive Smoothing,
Neuroimage 39(4), 1763-1773 (2008).
}
\details{
\tabular{ll}{
Package: \tab dti\cr
Version: \tab 0.5-7\cr
Date: \tab 2008-07-24\cr
Depends: \tab R (>= 2.5.0), adimpro, fmri, rgl\cr
License: \tab GPL (>=2)\cr
Copyright: \tab This package is
Copyright (C) 2005-2008 Weierstrass Institute for
Applied Analysis and Stochastics.\cr
URL: \tab \url{http://www.wias-berlin.de/projects/matheon_a3}\cr
Packaged: \tab Mon Oct 13 15:47:57 2008; tabelow\cr
Built: \tab R 2.7.1; i686-pc-linux-gnu; 2008-10-13 15:49:15; unix\cr
}

Index:
\preformatted{
clipdti-methods         Methods for Function 'clipdti' in Package 'dti'
dti-class               Class "dti"
dti.smooth              Structural Adaptive Smoothing
dtiData                 Read Diffusion Weighted Data from Image File
dtiData-class           Class "dtiData"
dtiIndices-class        Class "dtiIndices"
dtiIndices-methods      Methods for Function 'dtiIndices' in Package
                        'dti'
dtiTensor-class         Class "dtiTensor"
dtiTensor-methods       Methods for Function 'dtiTensor' in Package
                        'dti'
extract-methods         Methods for Function 'extract' in Package 'dti'
medinria2tensor         Read Tensor Data from NIFTI File
plot-methods            Methods for Function 'plot' in Package 'dti'
print-methods           Methods for Function 'print' in Package 'dti'
readDWIdata             Read Diffusion Weighted Data from Directory
show3d-methods          Methods for Function 'show3d' in Package 'dti'
summary-methods         Methods for Function 'summary' in Package 'dti'
tensor2medinria         Write Tensor Data as NIFTI File
}
}
\author{
Karsten Tabelow <tabelow@wias-berlin.de>
Joerg Polzehl <polzehl@wias-berlin.de>

Maintainer: Karsten Tabelow <tabelow@wias-berlin.de>
}
\keyword{ package }
\examples{
  \dontrun{demo(dti_art)}
}
