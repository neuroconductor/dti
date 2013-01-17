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

setMethod("dkiTensor", "dtiData",
          function(object, 
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

            ## define the design matrix of the estimation problem
            ind <- object@bvalue != 0
            xxx <- dkiDesign( object@gradient[ , ind])
            ## WHAT HAPPENS, IF BVALUE IS NOT ZERO, CAN WE EXTENT THE DESIGN MATRIX?
            
            ## Tabesh Eq. [11, 14, 15]
            Tabesh_AD <- xxx[ , 1:6]
            Tabesh_AK <- xxx[ , 7:21]
            
            ## Tabesh Eq. [10]
            Tabesh_A <- cbind(sweep( Tabesh_AD, 1, - object@bvalue[ ind], "*"),
                              sweep( Tabesh_AK, 1, object@bvalue[ ind]^2/6, "*"))

            ## Tabesh Eq. [13]
            Tabesh_C <- rbind(cbind( Tabesh_AD, matrix(0, 200, 15)), cbind( matrix(0, 200, 6), Tabesh_AK), cbind( 3/max(object@bvalue)*Tabesh_AD, - Tabesh_AK))

            ## container for estimation results
            ddim <- object@ddim  
            D <- array( 0, dim = c( 21, ddim))
            
            if ( method == "CLLS-QP") {

              if ( !require(quadprog)) return( "dkiTensor: did not find package quadprog, please install for the CLLS-QP method")

              ## scaling factor to make it numerically stable
              ## choice is arbitrary
              E <-  1e-3
              Dmat <- E^2 * t( Tabesh_A) %*% Tabesh_A
              Amat <- t( Tabesh_C)

              ## go through the dataset
              if (verbose) pb <- txtProgressBar(0, ddim[1], style = 3)
              for ( i in 1:ddim[1]) {
                if (verbose) setTxtProgressBar(pb, i)
                for ( j in 1:ddim[2]) {
                  for ( k in 1:ddim[3]) {
                    ## cat( "voxel", i, j, k, "\n")
				  
                    ## Tabesh Eq. [12]
                    Tabesh_B <- log( dwiobj@si[ i, j, k, -dwiobj@s0ind] / mean( dwiobj@si[ i, j, k, dwiobj@s0ind]))
                    ## PERFORM SOME TEST BEFORE IT!!
                    
                    ## QP solution
                    dvec <- E * as.vector( t( Tabesh_A) %*% Tabesh_B)
                    resQP <- solve.QP( Dmat, dvec, Amat)
                    D[ , i, j, k] <- E * resQP$solution[ c( 1, 4, 5, 2, 6, 3, 7:21)] # re-order DT estimate to comply with dtiTensor
                  }
                }
              }
              if (verbose) close(pb)
            }

            if ( verbose) cat( "dkiTensor: finished estimation", format( Sys.time()), "\n")	 

            dim(D) <- c( 21, prod(ddim))
            Dapp <- Tabesh_AD %*% D[ c( 1, 4, 6, 2, 3, 5), ]
            dim(D) <- c( 21, ddim)
            
            ## CHECK THESE!
            th0 <- apply( dwiobj@si, 1:3, mean)
            sigma2 <- array( 0, dim = ddim)
            scorr <- array( 0, dim =c(3,3,3))
            bw <- c( 0, 0, 0)
            index <- 1
            scale <- 1
            mask <- th0 > object@level
            
            if ( verbose) cat( "dkiTensor: exiting function", format( Sys.time()), "\n")	 
	  
	  
            invisible( new("dkiTensor",
                           call  = args,
                           D     = D,
                           th0   = th0,
                           sigma = sigma2,
                           scorr = scorr, 
                           bw = bw, 
                           mask = mask,
                           hmax = 1,
                           gradient = object@gradient,
                           bvalue = object@bvalue,
                           btb   = xxx,
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

setMethod("dkiIndices", "dkiTensor",
          function(object, 
                   mc.cores = getOption( "mc.cores", 2L),
                   verbose = FALSE) {
			 
            if ( verbose) cat( "dkiTensor: entering function", format( Sys.time()), "\n")
            
            ## call history  
            args <- c( object@call, sys.call(-1))

            ## we need this for all the arrays
            ddim <- object@ddim
            nvox <- prod( ddim)
            nvox0 <- sum( object@mask)

            ## perform the DTI indices calculations
            D <- object@D
            z <- dtiind3D( object@D[ 1:6, , , ], object@mask, mc.cores = mc.cores, verbose = verbose)

            if ( verbose) cat( "dkiTensor: DTI indices calculated", format( Sys.time()), "\n")

            ## perform the DKI indices calculations
            ## DO WE NEED THIS?
            if ( mc.cores > 1) {
              mc.cores.old <- setCores( , reprt = FALSE)
              setCores( mc.cores)
            }
            
            dim( D) <- c( 21, nvox)
            andir <- matrix( 0, 9, nvox)
            lambda <- matrix( 1e20, 3, nvox)
            
            t1 <- Sys.time()

            zz <- .Fortran("dti3DevAll",
                           as.double( D[ 1:6, object@mask]),
                           as.integer( nvox0),
                           andir = double( 9*nvox0),
                           evalues = double( 3*nvox0),
                           DUP = FALSE,
                           PACKAGE = "dti")[ c( "andir", "evalues")]
            andir[ , object@mask] <- zz$andir
            lambda[ , object@mask] <- zz$evalues 

            t2 <- Sys.time()
            if ( verbose) cat( "dkiTensor: calculation took ", difftime( t2, t1), attr(difftime( t2, t1), "units"), " for", nvox0, "voxel\n")
           
            if ( mc.cores > 1) setCores( mc.cores.old, reprt = FALSE)
            dim( andir) <- c( 3, 3, nvox)
            
            ind <- object@bvalue != 0
            xxx <- dkiDesign( object@gradient[ , ind])
            Tabesh_AD <- xxx[ , 1:6]
            Tabesh_AK <- xxx[ , 7:21]
            ## Note: the DT entries in D are re-ordered compared to Tabesh_AD for backward comp.
            ## Maybe we dont need this!
            Dapp <- Tabesh_AD %*% D[ c( 1, 4, 6, 2, 3, 5), ]
            Kapp <- (Tabesh_AK %*% D[ 7:21, ]) / Dapp^2
            ## remove pathological values!
            Kapp[ Kapp < 0] <- 0
            Kapp[ Dapp <= 0] <- 0
            ## define mean kurtosis
            mk <- apply( Kapp, 2, mean)
            
            x1 <- dkiDesign( andir[ , 1, ])
            x2 <- dkiDesign( andir[ , 2, ])
            x3 <- dkiDesign( andir[ , 3, ])

            ## cannot allocate memory for the following:
            ## w1111 <- diag( x1[ , 7:21] %*% D[ 7:21, ])
            w1111 <- numeric( nvox)
            for ( i in 1:nvox) w1111[ i] <- x1[ i, 7:21] %*% D[ 7:21, i]
            w2222 <- numeric( nvox)
            for ( i in 1:nvox) w2222[ i] <- x2[ i, 7:21] %*% D[ 7:21, i]
            w3333 <- numeric( nvox)
            for ( i in 1:nvox) w3333[ i] <- x3[ i, 7:21] %*% D[ 7:21, i]
            
            k1 <- w1111 / lambda[ 1, ]^2
            k2 <- w2222 / lambda[ 2, ]^2
            k3 <- w3333 / lambda[ 3, ]^2

            kbar <- ( k1 + k2 + k3) / 3
            kaxial <- k1
            kradial <- ( k2 + k3) / 2
            fak <- sqrt( 3/2 *( (k1-kbar)^2 + (k2-kbar)^2 + (k3-kbar)^2) / ( k1^2 + k2^2 + k3^2))

            ## we finally got k1, k2, k3, kaxial, kradial, mk, fak
            dim( k1) <- ddim
            dim( k2) <- ddim
            dim( k3) <- ddim
            dim( mk) <- ddim
            dim( kaxial) <- ddim
            dim( kradial) <- ddim
            dim( fak) <- ddim
            
            if ( verbose) cat( "dkiTensor: DKI indices calculated", format( Sys.time()), "\n")
                      
            if ( verbose) cat( "dkiTensor: exiting function", format( Sys.time()), "\n")

            invisible( new("dkiIndices",
                           call = args,
                           fa = array( z$fa, ddim),
                           ga = array( z$ga, ddim),
                           md = array( z$md, ddim),
                           andir = array( z$andir, c( 3, ddim)),
                           bary = array( z$bary, c( 3, ddim)),
                           k1 = k1,
                           k2 = k2,
                           k3 = k3,
                           mk = mk,
                           kaxial = kaxial,
                           kradial = kradial,
                           fak = fak,
                           gradient = object@gradient,
                           bvalue = object@bvalue,
                           btb   = object@btb,
                           ngrad = object@ngrad,
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
                           source = object@source)
                      )
})

dkiDesign <- function( gradients) {
	
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

