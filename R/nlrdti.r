testnlrdti <- function(s,gradient,th0=0,D=c(1,0,0,1,0,1),niter=10,method=1){
nb <- length(s)
if (dim(gradient)[2]==3) gradient <- t(gradient)
btb <- create.designmatrix.dti(gradient, bvalue=1)
if(method==1){
z <- .Fortran("solvedti",
              as.integer(s),
              as.integer(nb),
              as.double(btb),
              th0=as.double(th0),
              D=as.double(D),
              Varth=double(28),
              res=double(nb),
              as.integer(niter),
              as.double(1e-4),
              double(nb),
              rss=double(1),
              PACKAGE="dti",DUP=FALSE)[c("th0","D","Varth","res","rss")]
} else{
z <- .Fortran("solvedtb",
              as.integer(s),
              as.integer(nb),
              as.double(btb),
              th0=as.double(th0),
              D=as.double(D),
              Varth=double(28),
              res=double(nb),
              as.integer(niter),
              as.double(1e-4),
              double(nb),
              rss=double(1),
              PACKAGE="dti",DUP=FALSE)[c("th0","D","Varth","res","rss")]
}
z
}

create.voxel.data <- function(th0,D,gradient,sigma){
nb <- dim(b)[2]
if (dim(gradient)[2]==3) gradient <- t(gradient)
btb <- create.designmatrix.dti(gradient, bvalue=1)
si <- exp(-D%*%btb)*th0
as.integer(pmax(0,rnorm(si,si,pmin(th0/2.5,sigma))))
}
