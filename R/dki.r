###################################################################
#                                                                 #
# Diffusion Kurtosis Imaging (DKI)                                #
#                                                                 # 
# Literature: Jensen, J. H. et al. (2005), MRM 53, 1432-1440      #
#             Jensen, J. H. et al. (2010), NMR Biomed 23, 698-710 #
#             Tabesh, A. et al. (2011), MRM 65, 823-836           #
#             Hui et al. (2008), NI 42, 122-134                   #
#                                                                 #
###################################################################

dkiTensor <- function(object,  ...) cat( "No DKI tensor calculation defined for this class:", class( object), "\n")

setGeneric("dkiTensor", function( object, ...) standardGeneric( "dkiTensor"))

setMethod( "dkiTensor", "dtiData",
function( object, 
          method = c( "ULLS", "CLLS-QP", "CLLS-H") , 
          mc.cores = getOption( "mc.cores", 1L),
          verbose = FALSE) {

	  if ( verbose) cat( "dkiTensor: entering function", format( Sys.time()), "\n")
	  
	  ## check method
	  method <- match.arg( method)
	  
	  ## call history
	  args <- c( object@call, sys.call(-1))
	  
	  ## check available number of cores for parallel computation
	  if( is.null( mc.cores)) mc.cores <- 1
	  mc.cores <- min( mc.cores, detectCores())
	  
	  if ( verbose) cat( "dkiTensor: exiting function", format( Sys.time()), "\n")	 
	  
	  ind <- bvalues != 0
	  xxx.alt <- dkiDesign.alt( gradients[ , ind])
	  
	  ## Tabesh Eq. [10]
	  Tabesh_A <- cbind( sweep( xxx.alt[ , 1:6], 1, - bvalues[ ind], "*"),
			  sweep( xxx.alt[ , 7:21], 1, bvalues[ ind]^2/6, "*"))
	  ## then we have
	  ## Tabesh Eq. [11, 14, 15]
	  Tabesh_AD <- xxx.alt[ , 1:6]
	  Tabesh_AK <- xxx.alt[ , 7:21]
	  
	  ## Tabesh Eq. [12]
	  Tabesh_B <- log(dwiobj@si[ 41, 20, 15, -dwiobj@s0ind] / mean( dwiobj@si[ 41, 20, 15, dwiobj@s0ind]))
	  
	  ## Tabesh Eq. [13]
	  Tabesh_C <- rbind(cbind( Tabesh_AD, matrix(0, 200, 15)), cbind( matrix(0, 200, 6), Tabesh_AK), cbind( 3/max(bvalues)*Tabesh_AD, - Tabesh_AK))
	  
	  ## Tabesh Eq. [16] with Kmin = 0
	  Tabesh_d <- rep( 0, 600)
	  
	  ## QP
	  require(quadprog)
	  E <- 1e-3
	  Dmat <- E^2 * t( Tabesh_A) %*% Tabesh_A
	  dvec <- E * as.vector( t( Tabesh_A) %*% Tabesh_B)
	  Amat <- t( Tabesh_C)
	  resQP <- solve.QP( Dmat, dvec, Amat)
	  E * resQP$solution 
	  dtiobj@D[ c( 1, 4, 6, 2, 3, 5), 41, 20, 15]
	  
	  ddim <- dwiobj@ddim
	  dkiobj <- dtiobj
	  
	  for ( i in 21:100) {
		  for ( j in 1:48) {
			  for ( k in 15) {
				  cat( "voxel", i, j, k, "\n")
				  
				  ## Tabesh Eq. [12]
				  Tabesh_B <- log(dwiobj@si[ i, j, k, -dwiobj@s0ind] / mean( dwiobj@si[ i, j, k, dwiobj@s0ind]))
				  
				  ## QP solution
				  dvec <- E * as.vector( t( Tabesh_A) %*% Tabesh_B)
				  resQP <- solve.QP( Dmat, dvec, Amat)
				  dkiobj@D[ , i, j, k] <- resQP$solution[ c( 1, 4, 5, 2, 6, 3)]
			  }
		  }
	  }
	  
	  
	  invisible( new("dkiTensor",
                         call  = args,
                         D     = D,
                         th0   = th0,
                         sigma = sigma2,
                         scorr = scorr$scorr, 
                         bw = scorr$bw, 
                         mask = mask,
                         hmax = 1,
                         gradient = object@gradient,
                         bvalue = object@bvalue,
                         btb   = object@btb,
                         ngrad = object@ngrad,
                         s0ind = object@s0ind,
                         replind = object@replind,
                         ddim  = object@ddim,
                         ddim0 = object@ddim0,
                         xind  = object@xind,
                         yind  = object@yind,
                         zind  = object@zind,
                         voxelext = object@voxelext,
                         level = object@level,
                         orientation = object@orientation,
                         rotation = object@rotation,
                         source = object@source,
                         outlier = index,
                         scale = scale,
                         method = method)
                    )
})


dkiIndices <- function( object, ...) cat( "No DKI indices calculation defined for this class:", class( object), "\n")

setGeneric( "dkiIndices", function( object, ...) standardGeneric( "dkiIndices"))

setMethod( "dtiIndices", "dtiTensor",
function( object, 
          mc.cores = getOption( "mc.cores", 2L)
          verbose = FALSE) {
			 
          if ( verbose) cat( "dkiTensor: entering function", format( Sys.time()), "\n")
			  
          ## call history  
	      args <- c( object@call, sys.call(-1))

          invisible( new("dkiIndices",
                         call = args,
                         fa = array(z$fa,ddim),
                         ga = array(z$ga,ddim),
                         md = array(z$md,ddim),
                         andir = array(z$andir,c(3,ddim)),
                         bary = array(z$bary,c(3,ddim)),
                         gradient = object@gradient,
                         bvalue = object@bvalue,
                         btb   = object@btb,
                         ngrad = object@ngrad, # = dim(btb)[2]
                         s0ind = object@s0ind,
                         ddim  = ddim,
                         ddim0 = object@ddim0,
                         voxelext = object@voxelext,
                         orientation = object@orientation,
                         rotation = object@rotation,
                         xind  = object@xind,
                         yind  = object@yind,
                         zind  = object@zind,
                         method = object@method,
                         level = object@level,
                         source= object@source)
                    )
})

dkiDesign.alt <- function( gradients) {
	
	## start with some checks
	dgrad <- dim( gradients)
	
	## do we have a matrix (check for class is not correct here, as array or dataframe would also do)
	if ( length(dgrad) != 2) stop( "dkiDesign: dimensionality of gradients is not 2!")
	
	## we expect the columns of gradients to be the gradient vectors
	if ( dgrad[1] != 3) stop( "dkiDesign: not a valid gradient matrix, expect gradient vectors as columns!")
	
	## now we have: - dgrad[2] is the number of gradients
	
	X <- matrix( 0, dgrad[2], 21)
	
	X[ , 1] <- gradients[ 1, ]^2 # D_11
	X[ , 2] <- gradients[ 2, ]^2 # D_22
	X[ , 3] <- gradients[ 3, ]^2 # D_33
	
	X[ , 4] <- 2 * gradients[ 1, ] * gradients[ 2, ] # D_12
	X[ , 5] <- 2 * gradients[ 1, ] * gradients[ 3, ] # D_13
	X[ , 6] <- 2 * gradients[ 2, ] * gradients[ 3, ] # D_23
	
	X[ , 7] <- gradients[ 1, ]^4 # V_1111
	X[ , 8] <- gradients[ 2, ]^4 # V_2222
	X[ , 9] <- gradients[ 3, ]^4 # V_3333
	
	X[ , 10] <- 4 * gradients[ 1, ]^3 * gradients[ 2, ] # V_1112
	X[ , 11] <- 4 * gradients[ 1, ]^3 * gradients[ 3, ] # V_1113
	X[ , 12] <- 4 * gradients[ 2, ]^3 * gradients[ 1, ] # V_2221
	X[ , 13] <- 4 * gradients[ 2, ]^3 * gradients[ 3, ] # V_2223
	X[ , 14] <- 4 * gradients[ 3, ]^3 * gradients[ 1, ] # V_3331
	X[ , 15] <- 4 * gradients[ 3, ]^3 * gradients[ 2, ] # V_3332
	
	X[ , 16] <- 6 * gradients[ 1, ]^2 * gradients[ 2, ]^2 # V_1122
	X[ , 17] <- 6 * gradients[ 1, ]^2 * gradients[ 3, ]^2 # V_1133
	X[ , 18] <- 6 * gradients[ 2, ]^2 * gradients[ 3, ]^2 # V_2233
	
	X[ , 19] <- 12 * gradients[ 1, ]^2 * gradients[ 2, ] * gradients[ 3, ] # V_1123
	X[ , 20] <- 12 * gradients[ 2, ]^2 * gradients[ 1, ] * gradients[ 3, ] # V_2213
	X[ , 21] <- 12 * gradients[ 3, ]^2 * gradients[ 1, ] * gradients[ 2, ] # V_3312
	
	X
}

