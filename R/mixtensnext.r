pmixtn0b <- function(x,ngrad,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
  nvox <- length(x)/(ngrad+2+2*maxcomp)
  dim(x) <- c(ngrad+2+2*maxcomp,nvox)
  z <- .C(C_mixtrl0b,
          as.integer(nvox),#n1
          as.integer(x[(ngrad+2):(ngrad+1+maxcomp),]),#siind
          as.integer(x[-(1:(ngrad+1+maxcomp)),]),#wi
          as.integer(ngrad),#ngrad
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          mix     = double(maxcomp*nvox),#mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","mix")]
  lev <- matrix(0,2,nvox)
  lev[1,] <- (alpha+1)*lambda
  lev[2,] <- lambda
  rbind(z$order,z$sigma2,matrix(lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}

pmixtn1b <- function(x,ngrad,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
  nvox <- length(x)/(ngrad+2+2*maxcomp)
  dim(x) <- c(ngrad+2+2*maxcomp,nvox)
  z <- .C(C_mixtrl1b,
          as.integer(nvox),#n1
          as.integer(x[(ngrad+2):(ngrad+1+maxcomp),]),#siind
          as.integer(x[-(1:(ngrad+1+maxcomp)),]),#wi
          as.integer(ngrad),#ngrad
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          lambda  = double(nvox),#lambda_ret lambda_2
          mix     = double(maxcomp*nvox),#mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","lambda","mix")]
  lev <- matrix(0,2,nvox)
  lev[1,] <- (alpha+1)*z$lambda
  lev[2,] <- z$lambda
  rbind(z$order,z$sigma2,matrix(lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}

pmixtn2b <- function(x,ngrad,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
  nvox <- length(x)/(ngrad+2+2*maxcomp)
  dim(x) <- c(ngrad+2+2*maxcomp,nvox)
  z <- .C(C_mixtrl2b,
          as.integer(nvox),#n1
          as.integer(x[(ngrad+2):(ngrad+1+maxcomp),]),#siind
          as.integer(x[-(1:(ngrad+1+maxcomp)),]),#wi
          as.integer(ngrad),#ngrad
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          alpha   = double(nvox),#alpha_ret alpha=(lambda_1-lambda_2)/lambda_2
          lambda  = double(nvox),#lambda_ret lambda_2
          mix     = double(maxcomp*nvox),#mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")]
  lev <- matrix(0,2,nvox)
  lev[1,] <- (z$alpha+1)*z$lambda
  lev[2,] <- z$lambda
  rbind(z$order,z$sigma2,lev,matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}
pgetsiindbv <- function(x,grad,bv,nvico,dgrad,dgradi,isample,alpha,lambda,
                       maxcomp,maxc,nguess){
  #x contains
  # si:      x[1:nsi,]
  # sigma2:  x[nsi+1,]
  nvox <- dim(x)[2]
  nsi <- dim(x)[1]-1
  siind <- .Fortran(C_getsiibv,
                    as.double(x[1:nsi,]),
                    as.integer(nsi),
                    as.integer(nvox),
                    as.integer(maxcomp),
                    as.double(dgrad),
                    as.double(bv),
                    as.integer(nvico),
                    as.double(alpha),
                    as.double(lambda),
                    double(nsi*nvico),
                    as.integer(isample),
                    as.integer(nguess),
                    double(nsi),
                    double(nsi),#z0
                    double(nsi*(maxcomp+1)),
                    siind=integer((maxcomp+1)*nvox),
                    wi=double((maxcomp+1)*nvox),
                    krit=double(nvox),
                    as.integer(maxcomp+1),
                    PACKAGE="dti")[c("siind","krit","wi")]
  z <- matrix(0,2*maxcomp+3,nvox)
  z[2:(maxcomp+2),] <- siind$siind  ## vertex indices
  z[-(1:(maxcomp+2)),] <- siind$wi  ## weights
  z[1,] <- siind$krit  ## kriterion
  z
}
