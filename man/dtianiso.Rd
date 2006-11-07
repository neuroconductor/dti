\name{dtianiso}
\alias{dtianiso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Adaptive (anisotropic) smoothing of 3D diffusion tensor images. }
\description{ The function performs adaptive (anisotropic) smoothing of 3D diffusion tensor images.
}
\usage{
dtianiso(y, hmax=5, lambda=25, rho=1, graph = FALSE, slice = NULL, bvec = NULL, sigma2 = NULL, scorr = c(0.5, 0.5), mask = NULL, quant = 0.8, minanindex = NULL, zext = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ array of dimension \code{c(6,nx,ny,nz)} containing the estimated tensors in each voxel. 
           tensor components are stored as \code{c(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz)}.
  \item{hmax}{\code{hmax} Maximal bandwidth to be used. }
  \item{lambda}{\code{lambda} Scale parameter in the statistical penalty,  }
  \item{rho}{\code{rho} Regularization parameter for anisotropic neighborhoods. Large values lead
      to isotropic smoothing, for small values ellipsoids defined by the diffusion tensor are used. 
      The effect of regularization vanishes with reduced variability of the estimates.  }
  \item{graph}{\code{graph} Logical, if \code{graph==TRUE} anisotropy information and sum of weights
     are illustrated for the slice specified in \code{slice}. The behaviour is influenced by 
     \code{quant} and \code{minanindex} .}
  \item{slice}{ If\code{graph==TRUE} \code{slice} specifies the slice to show.   }
  \item{bvec}{ \code{bvec} matrix of gradient vectors applied when recording the DT image.
             Dimension is \code{c(3,ngradients)}. }
  \item{sigma2}{ \code{sigma2} contains an array of voxelwise variances for \code{y} }
  \item{scorr}{ \code{scorr} specifies the spatial correlation within slides. The correlation between slides
    is assumed to be zero. }
  \item{mask}{ \code{mask} Logical array of size \code{dim(y)[-1]}. Computations are restricted to
    voxel where \code{mask} is \code{TRUE}  }
  \item{quant}{ If \code{graph} and \code{is.null(minanindex)} then \code{quant} 
       determines the lower bound for the anisotropy index as the \code{quant}-quantile of estimated
       anisotropy indices.  }
  \item{minanindex}{If \code{graph}  \code{minanindex} allows to determine the lower bound for the anisotropy index. 
   only results at voxel with anisotropy index larger than this lower bound will be displayed. }
  \item{zext}{ \code{zext} specifies a factor by which the vertical extension (distance between slices) of a voxel 
   differs for the voxelsize within a slice. }
}
\details{
}
\value{ A list with components
theta=z$theta,bi=z$bi,anindex=z$anindex,andirection=z$andirection,mask=z$mask,call=args  
  If it is a LIST, use
  \item{theta }{The estimated (smoothed) diffusion tensor. Same dimension as \code{y}.}
  \item{bi }{ \code{bi[i,j,k]*InvCov} is an estimate of the inverse of the covariance matrix of theta in
         voxel \code{c(i,j,k)}. }
  \item{anindex }{Array of anisotropy indices. Dimension\code{dim(y)[-1]}. }
  \item{andirection }{Array containing the first eigenvectors of the estimated diffusion tensor. Dimension\code{c(3,dim(y)[-1])}.}
  \item{mask }{Logical mask used to restrict the computations to voxel wit \code{mask==TRUE}.}
  \item{InvCov }{\code{bi[i,j,k]*InvCov} is an estimate of the inverse of the covariance matrix of theta in
         voxel \code{c(i,j,k)}, \code{InvCov/sigma2[i,j,k]} the estimated inverse covariance matrix of the observed tensor.}
  \item{args }{Arguments of the call to \code{dtianiso}.}
  The returned object has an class attribute \code{"dti"}.
}
\references{   Polzehl and Spokoiny (2006). Propagation-Separation Approach for Local
  Likelihood Estimation. Probability Theory and Related Fields. 3 (135)
  335 - 362.

  Polzehl, J. and Tabelow, K. (2006). Adaptive smoothing of DT images, manuscript. 
}
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de }
\note{ }
\seealso{ \code{\link{andir.image}}}
\examples{
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
