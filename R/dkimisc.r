rotateKurtosis <- function( evec, KT, i = 1, j = 1, k = 1, l = 1) {
  
  if ( ( i < 1) | ( i > 3)) stop( "rotateKurtosis: index i out of range")
  if ( ( j < 1) | ( j > 3)) stop( "rotateKurtosis: index j out of range")
  if ( ( k < 1) | ( k > 3)) stop( "rotateKurtosis: index k out of range")  
  if ( ( l < 1) | ( l > 3)) stop( "rotateKurtosis: index l out of range")

  if ( length( dim( evec)) != 3) stop( "rotateKurtosis: dimension of direction array is not 3")
  if ( dim( evec)[1] != 3) stop( "rotateKurtosis: length of direction vector is not 3")
  if ( dim( evec)[2] != 3) stop( "rotateKurtosis: number of direction vectors is not 3")

  if ( length( dim( KT)) != 2) stop( "rotateKurtosis: dimension of kurtosis array is not 2")
  if ( dim( KT)[1] != 15) stop( "rotateKurtosis: kurtosis tensor does not have 15 elements")
  if ( dim( KT)[2] != dim( evec)[3]) stop( "rotateKurtosis: number of direction vectors does not match number of kurtosis tensors")
  
  nvox <- dim( KT)[2]
  
  ## we create the full symmetric kurtosis tensor for each voxel ...
  W <- defineKurtosisTensor( KT)
  ## ..., i.e., dim(W) == c( 3, 3, 3, 3, nvox) 
  
  Wtilde <- rep( 0, nvox)
  
  for ( ii in 1:3) {
    for ( jj in 1:3) {
      for ( kk in 1:3) {
        for ( ll in 1:3) {
          
          ## TODO: I am not sure about the order of i and ii, should it be evec[ i, ii, n]? Same for j, k, l
          Wtilde <- Wtilde + evec[ ii, i, ] * evec[ jj, j, ] * evec[ kk, k, ] * evec[ ll, l, ] * W[ ii, jj, kk, ll, ]   
#          Wtilde <- Wtilde + evec[ i, ii, ] * evec[ j, jj, ] * evec[ k, kk, ] * evec[ l, ll, ] * W[ ii, jj, kk, ll, ]   
          
        }
      }
    }
  }
  
  invisible( Wtilde)
  
}


defineKurtosisTensor <- function( DK) {
  
  W <- array( 0, dim = c( 3, 3, 3, 3, dim( DK)[2]))
  
  W[ 1, 1, 1, 1, ] <- DK[ 1, ]
  W[ 2, 2, 2, 2, ] <- DK[ 2, ]
  W[ 3, 3, 3, 3, ] <- DK[ 3, ]
  
  W[ 1, 1, 1, 2, ] <- DK[ 4, ]
  W[ 1, 1, 2, 1, ] <- DK[ 4, ]
  W[ 1, 2, 1, 1, ] <- DK[ 4, ]
  W[ 2, 1, 1, 1, ] <- DK[ 4, ]
  
  W[ 1, 1, 1, 3, ] <- DK[ 5, ]
  W[ 1, 1, 3, 1, ] <- DK[ 5, ]
  W[ 1, 3, 1, 1, ] <- DK[ 5, ]
  W[ 3, 1, 1, 1, ] <- DK[ 5, ]
  
  W[ 2, 2, 2, 1, ] <- DK[ 6, ]
  W[ 2, 2, 1, 2, ] <- DK[ 6, ]
  W[ 2, 1, 2, 2, ] <- DK[ 6, ]
  W[ 1, 2, 2, 2, ] <- DK[ 6, ]
  
  W[ 2, 2, 2, 3, ] <- DK[ 7, ]
  W[ 2, 2, 3, 2, ] <- DK[ 7, ]
  W[ 2, 3, 2, 2, ] <- DK[ 7, ]
  W[ 3, 2, 2, 2, ] <- DK[ 7, ]
  
  W[ 3, 3, 3, 1, ] <- DK[ 8, ]
  W[ 3, 3, 1, 3, ] <- DK[ 8, ]
  W[ 3, 1, 3, 3, ] <- DK[ 8, ]
  W[ 1, 3, 3, 3, ] <- DK[ 8, ]
  
  W[ 3, 3, 3, 2, ] <- DK[ 9, ]
  W[ 3, 3, 2, 3, ] <- DK[ 9, ]
  W[ 3, 2, 3, 3, ] <- DK[ 9, ]
  W[ 2, 3, 3, 3, ] <- DK[ 9, ]
  
  W[ 1, 1, 2, 2, ] <- DK[ 10, ]
  W[ 1, 2, 1, 2, ] <- DK[ 10, ]
  W[ 1, 2, 2, 1, ] <- DK[ 10, ]
  W[ 2, 1, 2, 1, ] <- DK[ 10, ]
  W[ 2, 1, 1, 2, ] <- DK[ 10, ]
  W[ 2, 2, 1, 1, ] <- DK[ 10, ]
  
  W[ 1, 1, 3, 3, ] <- DK[ 11, ]
  W[ 1, 3, 1, 3, ] <- DK[ 11, ]
  W[ 1, 3, 3, 1, ] <- DK[ 11, ]
  W[ 3, 1, 3, 1, ] <- DK[ 11, ]
  W[ 3, 1, 1, 3, ] <- DK[ 11, ]
  W[ 3, 3, 1, 1, ] <- DK[ 11, ]
  
  W[ 2, 2, 3, 3, ] <- DK[ 12, ]
  W[ 2, 3, 2, 3, ] <- DK[ 12, ]
  W[ 2, 3, 3, 2, ] <- DK[ 12, ]
  W[ 3, 2, 3, 2, ] <- DK[ 12, ]
  W[ 3, 2, 2, 3, ] <- DK[ 12, ]
  W[ 3, 3, 2, 2, ] <- DK[ 12, ]
  
  W[ 1, 1, 2, 3, ] <- DK[ 13, ]
  W[ 1, 1, 3, 2, ] <- DK[ 13, ]
  W[ 3, 1, 1, 2, ] <- DK[ 13, ]
  W[ 2, 1, 1, 3, ] <- DK[ 13, ]
  W[ 2, 3, 1, 1, ] <- DK[ 13, ]
  W[ 3, 2, 1, 1, ] <- DK[ 13, ]
  W[ 1, 3, 1, 2, ] <- DK[ 13, ]
  W[ 1, 2, 1, 3, ] <- DK[ 13, ]
  W[ 3, 1, 2, 1, ] <- DK[ 13, ]
  W[ 2, 1, 3, 1, ] <- DK[ 13, ]
  W[ 1, 2, 3, 1, ] <- DK[ 13, ]
  W[ 1, 3, 2, 1, ] <- DK[ 13, ]
  
  W[ 2, 2, 1, 3, ] <- DK[ 14, ]
  W[ 2, 2, 3, 1, ] <- DK[ 14, ]
  W[ 3, 2, 2, 1, ] <- DK[ 14, ]
  W[ 1, 2, 2, 3, ] <- DK[ 14, ]
  W[ 1, 3, 2, 2, ] <- DK[ 14, ]
  W[ 3, 1, 2, 2, ] <- DK[ 14, ]
  W[ 3, 2, 1, 2, ] <- DK[ 14, ]
  W[ 2, 3, 2, 1, ] <- DK[ 14, ]
  W[ 1, 2, 3, 2, ] <- DK[ 14, ]
  W[ 2, 1, 2, 3, ] <- DK[ 14, ]
  W[ 2, 1, 3, 2, ] <- DK[ 14, ]
  W[ 2, 3, 1, 2, ] <- DK[ 14, ]
  
  W[ 3, 3, 1, 2, ] <- DK[ 15, ]
  W[ 3, 3, 2, 1, ] <- DK[ 15, ]
  W[ 2, 3, 3, 1, ] <- DK[ 15, ]
  W[ 1, 3, 3, 2, ] <- DK[ 15, ]
  W[ 1, 2, 3, 3, ] <- DK[ 15, ]
  W[ 2, 1, 3, 3, ] <- DK[ 15, ]
  W[ 1, 3, 2, 3, ] <- DK[ 15, ]
  W[ 2, 3, 1, 3, ] <- DK[ 15, ]
  W[ 3, 2, 3, 1, ] <- DK[ 15, ]
  W[ 3, 1, 3, 2, ] <- DK[ 15, ]
  W[ 3, 1, 2, 3, ] <- DK[ 15, ]
  W[ 3, 2, 1, 3, ] <- DK[ 15, ]
  
  invisible( W)
}

## TODO: For efficiency we might combine both functions
##       Many doubled calculations
## A strange small discontinuity at singularities for F1 remains

kurtosisFunctionF1<- function( l1, l2, l3) {
  
  require( gsl)
  ## Tabesh et al. Eq. [28], [A9, A12]
  ## this function is defined without 9*MD^2!!
  
  ## consider removable singularities!! 
  F1 <- numeric( length( l1))
  
  ind12 <- abs( l1 - l2) > 1e-10 ## l1 != l2
  ind13 <- abs( l1 - l3) > 1e-10 ## l1 != l3
  
  ind <- ( ind12) & ( ind13)
  if ( any( ind)) F1[ ind] <- ( ellint_RF( l1[ ind]/l2[ ind], l1[ ind]/l3[ ind], 1) * sqrt( l2[ ind]*l3[ ind]) / l1[ ind] + ellint_RD( l1[ ind]/l2[ ind], l1[ ind]/l3[ ind], 1) * ( 3* l1[ ind]^2 - l1[ ind]*l2[ ind] - l1[ ind]*l3[ ind] - l2[ ind]*l3[ ind]) / (3*l1[ ind]*sqrt( l2[ ind]*l3[ ind])) - 1) / 2 / ( l1[ ind]-l2[ ind]) / ( l1[ ind]-l3[ ind])
  
  ind1 <- ( !ind12) & ( ind13)
  if ( any( ind1)) F1[ ind1] <- kurtosisFunctionF2( l2[ ind1], l1[ ind1], l1[ ind1]) / 2 
  
  ind2 <- ( ind12) & ( !ind13)
  if ( any( ind2)) F1[ ind2] <- kurtosisFunctionF2( l3[ ind2], l1[ ind2], l1[ ind2]) / 2 
  
  ind3 <- ( !ind12) & ( !ind13)
  if ( any( ind3)) F1[ ind3] <- 9/5/( l1[ ind3] + l2[ ind3] + l3[ ind3])^2
  
  F1
}

kurtosisFunctionF2 <- function( l1, l2, l3) {
  
  require( gsl)
  ## Tabesh et al. Eq. [28], [A10, A11, A12]
  ## this function is defined without 9 * MD^2!!
  
  alpha <- function(x) 1/sqrt(abs(x)) * atan(sqrt(abs(x)))
  
  ## consider removable singularities!!
  F2 <- numeric( length( l1))
  
  ind23 <- abs( l2 - l3) > 1e-10 ## l2 != l3
  ind12 <- abs( l1 - l2) > 1e-10 ## l1 != l2
  
  if ( any( ind23)) F2[ ind23] <- 3 * ( ellint_RF( l1[ ind23]/l2[ ind23], l1[ ind23]/l3[ ind23], 1) * (l2[ ind23]+l3[ ind23]) / sqrt( l2[ ind23]*l3[ ind23]) + ellint_RD( l1[ ind23]/l2[ ind23], l1[ ind23]/l3[ ind23], 1) * (2*l1[ ind23]-l2[ ind23]-l3[ ind23]) / 3/sqrt( l2[ ind23]*l3[ ind23]) - 2 ) / (l2[ ind23] - l3[ ind23]) / (l2[ ind23] - l3[ ind23])
  
  ind1 <- ( !ind23) & ( ind12)
  if ( any( ind1)) F2[ ind1] <- 54 * (l1[ ind1]+2*l3[ ind1])^2/144/l3[ ind1]^2/(l1[ ind1]-l3[ ind1])^2 *(l3[ ind1]*(l1[ ind1]+2*l3[ ind1])+ l1[ ind1]*(l1[ ind1]-4*l3[ ind1])*alpha(1-l1[ ind1]/l3[ ind1]))/( l1[ ind1] + l2[ ind1] + l3[ ind1])^2
  
  ind2 <- ( !ind23) & ( !ind12)
  if ( any( ind2)) F2[ ind2] <- 54/15/( l1[ ind2] + l2[ ind2] + l3[ ind2])^2
  
  F2
}

pseudoinverseSVD <- function( xxx) {
  svdresult <- svd( xxx)
  svdresult$v %*% diag( 1 / svdresult$d) %*% t( svdresult$u)
}
