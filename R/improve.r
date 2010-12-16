improve <- function(mixtensobj,dwiobj, maxcomp=3,  p=40, method="mixtensor", reltol=1e-6, maxit=5000,ngc=1000, optmethod="BFGS", nguess=100*maxcomp^2,msc="BIC",pen=1,where=NULL){
#
#  uses  S(g)/s_0 = w_0 exp(-l_1) +\sum_{i} w_i exp(-l_2-(l_1-l_2)(g^T d_i)^2)
#
#  
#
  set.seed(1)
  theta <- .5
  maxc <- .866
  args <- sys.call(-1)
  args <- c(mixtensobj@call,args)
  ngrad <- mixtensobj@ngrad
  ddim <- mixtensobj@ddim
  mask <- mixtensobj@mask
  if(is.null(where)||any(dim(where)!=ddim[1:3])) where <- mask
  where <- where & mask
  s0ind <- mixtensobj@s0ind
  vext <- mixtensobj@voxelext
  if(any(ddim != dwiobj@ddim)){
     warning("incompatible Mixtensor- and dwiData objects, returning original Mixtensor-object")
     return(mixtensobj)
  }
  ns0 <- length(s0ind)
  ngrad0 <- ngrad - ns0
  if(5*(1+3*maxcomp)>ngrad0){
#     maxcomp <- max(1,trunc((ngrad0-5)/15))
     cat("Maximal number of components reduced to", maxcomp,"due to insufficient
          number of gradient directions\n")
  }
#
#  First tensor estimates to generate eigenvalues and -vectors
#
  prta <- Sys.time()
  gc()
  z <- extract(mixtensobj,c("order","mix","andir"))
  nmix <- mix <- z$mix
  norder <- order <- z$order
  nandir <- andir <- z$andir
  maxorder <- dim(andir)[2]
  maxcomp <- min(maxcomp,maxorder) 
  rm(z)
  gc()
  nsigma2 <- sigma2 <- mixtensobj@sigma
  norient <- orient <- mixtensobj@orient
  nlev <- lev <-mixtensobj@ev
  mask <- mixtensobj@mask
#  cat("Start search outlier detection at",format(Sys.time()),"\n")
#
#  replace physically meaningless S_i by mena S_0 values
#
#  z <- .Fortran("outlier1",#misc.f 
#                as.double(dwiobj@si),
#                as.logical(where),
#                as.integer(prod(ddim)),
#                as.integer(ngrad),
#                as.logical((1:ngrad)%in%s0ind),
#                as.integer(ns0),
#                si=integer(prod(ddim)*ngrad),
#                index=integer(prod(ddim)),
#                lindex=integer(1),
#                DUP=FALSE,
#                PACKAGE="dti")[c("si","index","lindex")]
#  cat("End search outlier detection at",format(Sys.time()),"\n")
#  si <- array(as.integer(z$si),c(ddim,ngrad))
#  rm(z)
#  gc()
#  gc()
  cat("Start generating auxiliary objects",format(Sys.time()),"\n")
#
#  compute mean S_0, s_i/S_0 (siq), var(siq) and mask
#
  sii <- matrix(dwiobj@si,prod(ddim[1:3]),ngrad)[where,]      
  z <- .Fortran("sweepimp",# mixtens.f
                as.integer(sii[,-s0ind]),
                as.integer(sii[,s0ind]),
                as.integer(sum(where)),
                as.integer(ns0),
                as.integer(ngrad0),
                siq=double(sum(where)*ngrad0),
                s0=double(sum(where)),
                DUPL=FALSE,
                PACKAGE="dti")[c("siq","s0")]
  rm(sii)
  s0 <- array(0,ddim[1:3])
  s0[where] <- z$s0
  siq <- matrix(0,prod(ddim[1:3]),ngrad0)
  siq[where,] <- z$siq
  dim(siq) <- c(ddim[1:3],ngrad0)
  rm(z)
  gc()
  gc()
  npar <- if(method=="mixtensor") 1+3*(0:maxcomp) else c(1,2+3*(1:maxcomp))
#
#   compute penalty for model selection, default BIC
#
  penIC <- switch(msc,"AIC"=2*npar/ngrad0,"BIC"=log(ngrad0)*npar/ngrad0,
                  "AICC"=(1+npar/ngrad0)/(1-(npar+2)/ngrad0),
                  log(ngrad0)*npar/ngrad0)
  cat("End generating auxiliary objects",format(Sys.time()),"\n")
#
#  avoid situations where si's are larger than s0
#
  grad <- t(dwiobj@gradient[,-s0ind])
#
#   determine initial estimates for orientations 
#
  n1 <- ddim[1]
  n2 <- ddim[2]
  n3 <- ddim[3]
  igc <- 0
  ingc <- 0
  ind3 <- -1:1
  prt0 <- Sys.time()
#
#   loop over voxel in volume
#  exclude indices at the sides of the cube (just to save index operations,
#  usually they are outside the mask anyway
#
  for(i3 in 2:(n3-1)) for(i2 in 2:(n2-1)) for(i1 in 2:(n1-1)){ # begin loop
     ordi <- order[i1,i2,i3]
     krit <- log(sigma2[i1,i2,i3])+penIC[ordi+1]
     if(where[i1,i2,i3]&(ordi<maxcomp)){ # begin mask
#   only analyze voxel within mask
     mc0 <- maxcomp
     ord <- mc0+1
     z <- .Fortran("imprparb",
                     as.integer(maxcomp),
                     as.integer(maxorder),
                     as.integer(order[i1+ind3,i2+ind3,i3+ind3]),
                     as.double(mix[,i1+ind3,i2+ind3,i3+ind3]),
                     as.double(andir[,,i1+ind3,i2+ind3,i3+ind3]),
                     as.double(orient[,,i1+ind3,i2+ind3,i3+ind3]),
                     as.double(lev[,i1+ind3,i2+ind3,i3+ind3]),
                     as.double(vext),
                     param=numeric(2*maxcomp+1),
                     as.integer(2*maxcomp+1),
                     npar=integer(1),
                     DUPL=TRUE,
                     PACKAGE="dti")[c("param","npar")]   
     par <- z$param
#
#   these are the gradient vectors corresponding to minima in spherical coordinates
#
#
#  use AIC/ngrad0, BIC/ngrad0 or AICC/ngrad0 respectively
#
#        cat("i",i1,i2,i3,"par",par[1:z$npar],"npar",z$npar,"\n")
     mc0 <- (z$npar+1)/2

     for(k in mc0:1){ # begin model order
        if(k<ord) {
#
#  otherwise we would reanalyze a model
#
        if(method=="mixtensor"){
           lpar <- 2*k+1
#
#        cat("par",par[1:(2*k+1)],"pen",pen,"krit",krit,"\n")
           if(optmethod=="BFGS"){
                 z <- optim(par[1:(2*k+1)],mfunpl0,gmfunpl0,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method="BFGS",control=list(maxit=maxit,reltol=reltol))
           } else {
              z <- optim(par[1:(2*k+1)],mfunpl0,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method=optmethod,control=list(maxit=maxit,reltol=reltol))
           }
        } else if(method=="mixtensoriso"){
           lpar <- 2*k+1
#
           if(optmethod=="BFGS"){
                 z <- optim(par[1:(2*k+1)],mfunpl0,gmfunpl0,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method="BFGS",control=list(maxit=maxit,reltol=reltol))
                 z <- optim(z$par,mfunpl1,gmfunpl1,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method="BFGS",control=list(maxit=maxit,reltol=reltol))
           } else {
              z <- optim(par[1:(2*k+1)],mfunpl1,siq=siq[i1,i2,i3,],grad=grad,pen=pen,
                         method=optmethod,control=list(maxit=maxit,reltol=reltol))
           }
        }         
#        cat("opt-par",z$par,"value",z$value,"krit",krit,"\n")
# thats sum of squared residuals + penalties (w<0 or 0>th or or th > 8)
#
#   estimate of sigma from the best fitting model
#
        if(method=="mixtensor"){
            zz <- mfunplwghts0(z$par[1:lpar],siq[i1,i2,i3,],grad,pen)
        } else if (method=="mixtensoriso"){
            zz <- mfunplwghts1(z$par[1:lpar],siq[i1,i2,i3,],grad,pen)
        }
        value <- zz$value 
        ord <- zz$ord
#  replace sigmai by best variance estimate from currently best model
        if(any(zz$lev<0)||ord<k){
           ttt <- krit
#   parameters not interpretable reduce order
        } else {
           si2new <- value/(ngrad0-3*ord-1)
           if(si2new<1e-15){
               cat(i1,i2,i3,ord,si2new,"\n")
               si2new <- 1e-15
               }
           ttt <- log(si2new)+penIC[1+ord]
           par <- zz$par
        }
#        cat("par",par,"value",value,"ord",ord,"w",zz$mix,"lev",zz$lev,"\n")
#
#     use directions corresponding to largest weights as initial directions
#
        if(ttt < krit) {
           krit <- ttt
           norder[i1,i2,i3] <- ord
           nlev[,i1,i2,i3] <- zz$lev
           nmix[,i1,i2,i3] <- if(ord==maxorder) zz$mix else c(zz$mix,rep(0,maxorder-ord))
           norient[,1:ord,i1,i2,i3] <- zz$orient
           nsigma2[i1,i2,i3] <- si2new
       }
     }
   } # end model order#
    if(igc<ngc){
       igc <- igc+1
    } else {
       igc <- 1
       ingc <- ingc+1
       prt1 <- Sys.time()
       gc()
       cat("Nr. of voxel",ingc*ngc,"time elapsed:",format(difftime(prt1,prta),digits=3),"remaining time:",
            format(difftime(prt1,prt0)/(ingc*ngc)*(sum(where)-ingc*ngc),digits=3),"\n")
    }
  }# end mask
  }# end loop
  invisible(new("dwiMixtensor",
                model = "homogeneous_prolate",
                call   = args,
                ev     = nlev,
                mix    = nmix,
                orient = norient,
                order  = norder,
                p      = mixtensobj@p,
                th0    = mixtensobj@th0,
                sigma  = nsigma2,
                scorr  = mixtensobj@scorr, 
                bw     = mixtensobj@bw, 
                mask   = mixtensobj@mask,
                hmax   = mixtensobj@hmax,
                gradient = mixtensobj@gradient,
                btb    = mixtensobj@btb,
                ngrad  = mixtensobj@ngrad, # = dim(btb)[2]
                s0ind  = mixtensobj@s0ind,
                replind = mixtensobj@replind,
                ddim   = mixtensobj@ddim,
                ddim0  = mixtensobj@ddim0,
                xind   = mixtensobj@xind,
                yind   = mixtensobj@yind,
                zind   = mixtensobj@zind,
                voxelext = mixtensobj@voxelext,
                level = mixtensobj@level,
                orientation = mixtensobj@orientation,
                rotation = mixtensobj@rotation,
                source = mixtensobj@source,
                outlier = mixtensobj@outlier,
                scale = mixtensobj@scale,
                method = mixtensobj@method)
            )
   }

