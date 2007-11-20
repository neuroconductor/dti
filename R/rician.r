rician <- function(y,theta,sigma2,w,eps){
   n <- length(y)
   theta <- .Fortran("ricefix",
                 as.double(y),
                 theta=as.double(theta),
                 as.double(sigma2),
                 as.double(w),
                 as.integer(n),
                 as.double(eps),
                 DUP=TRUE,
                 PACKAGE="dti")$theta
   vtheta <- .Fortran("ricevar",
                 as.double(y),
                 as.double(theta),
                 as.double(sigma2),
                 as.double(w),
                 as.integer(n),
                 vtheta=double(1),
                 DUP=TRUE,
                 PACKAGE="dti")$vtheta
list(theta=theta,vtheta=vtheta)
}
