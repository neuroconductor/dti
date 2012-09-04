Qlm <- function(l1,l2,a,m1,m2,b){
    sqrt((2*l1+1)*(2*l2+1)*(2*a+1)/4/pi)*
         coupling_3j(2*l1,2*l2,2*a,0,0,0)*
         coupling_3j(2*l1,2*l2,2*a,2*m1,2*m2,2*b)
}
kappan <- function(n,zeta){
sqrt(2/zeta^1.5*gamma(1:(n+1))/gamma(0:n+3/2))
}
hnni <- function(n){
lni <- function(n,i) (-1)^i*gamma(n+1.5)/gamma(n-i+1)/gamma(i+1.5)/gamma(i+1)
hi <- array(0,c(n+1,n+1,2*n+1))
for(n1 in 0:n){ 
   ln1i <- lni(n1,0:n1)
   for(n2 in 0:n){
      ln2i <- lni(n2,0:n2)
      for(i in 0:n1) for(j in 0:n2){
         k <- i+j
         hi[n1+1,n2+1,k+1] <- hi[n1+1,n2+1,k+1] + ln1i[i+1]*ln2i[j+1]
         }
      }
   }
hi
}
g1f1 <- function(n,L,zeta,q=1){
g1 <- array(0,c(L+1,2*n+1))
pi2zetaq <- pi^2*zeta*q^2
for(i in 0:(2*n)) for(j in (0:L)) {
   g1[j+1,i+1] <- gamma(j+i+1.5)*hyperg_1F1(i+j+1.5, 2*j+1.5, - pi2zetaq)
}
g1
}
Inna <- function(n,L,zeta,q=1){
Inn <- array(0,c(L+1,n+1,n+1))
hi <- hnni(n)
g1 <- g1f1(n,L,zeta,q=1)
kn <- kappan(n,zeta)
for(a in 0:L) {
   za <- zeta^(0:L+1.5)*pi^(2*L+.5)*q^(2*(0:L))/4/gamma(2*(0:L)+1.5)
   for(i in 0:n) for(j in 0:n){
       for(k in 0:(i+j)) {
          Inn[,i+1,j+1] <- Inn[,i+1,j+1]+ za*hi[i+1,j+1,k+1]*g1[,k+1] 
          }
       }
   }
Inn
}
Yab <- function(order,gradients){
  order <- as.integer(max(0,order))
  if(order%%2==1){
    warning("maximum order needs to be even, increase order by one")
    order <- order+1
  } 

  # calculate spherical angles theta and phi corresponding to the gradients
  n <- dim(gradients)[2]
  theta <- phi <- numeric(n)
  for( i in 1:n){
    angles <- sphcoord(gradients[,i])
    theta[i] <- angles[1]
    phi[i] <-  angles[2]
  }

  # values of SH on specified spherical angles 
  getsphericalharmonicseven(order,theta,phi)# first index spharmon (a,b)-index 
  #second index - gradient
}
Qlmmat <- function(order){
# l in (0:(L%/%2))*2 
order <- as.integer(max(0,order))
if(order%%2==1){
warning("maximum order needs to be even, increase order by one")
order <- order+1
} 
kseq <- seq(0,order,2)
nind <- (order+1)*(order+2)/2
kseq2 <- seq(0,2*order,2)
nind2 <- (2*order+1)*(order+1)
Qmat <- array(0,c(nind,nind,nind2))
i1 <- i2 <- i3 <- 1
for(l1 in kseq){
   m1seq <- seq(-l1,l1,1)
   for(m1 in m1seq){
      i2 <- 1
      for(l2 in kseq){
         m2seq <- seq(-l2,l2,1)      
         for(m2 in m2seq){
            i3 <- 1
            for(a in kseq2){
               bseq <- seq(-a,a,1)
               for(b in bseq){
                  if(m1+m2+b==0) Qmat[i1,i2,i3] <- Qlm(l1,l2,a,m1,m2,b)
#  otherwise the result is zero
                  i3 <- i3+1
                  }
               }               
               i2 <- i2+1
            }
         }
         i1 <- i1+1
      }
   }
Qmat
}

Kofgrad <- function(grad,zeta,N=1,L=2,q=1){
#
#   compute the kernel K(q|\zeta) for all gradient directions
#   according to Jian Cheng "Nonnegative Definite EAP and ODF Estimation ...
# 
require(gsl)
ngrad <- dim(grad)[2]
Np1 <- N+1
Lind <- (L+1)*(L+2)/2
Lind2 <- (2*L+1)*(L+1)
a <- NULL
for(i in 0:(L)) a <- c(a,rep(2*i,4*i+1))
Kofg <- array(0,c(Lind,Lind,N+1,N+1,ngrad))
Yabm <- Yab(2*L,grad)# dim Lind2, ngrad
Qlmm <- Qlmmat(L)# dim Lind,Lind,Lind2
Inn <- Inna(N,L,zeta,q=1)# dim L+1 (even values only) , N+1, N+1
for(i1 in 1:Lind) for(i2 in 1:i1)
for(n1 in 1:Np1) for(n2 in 1:n1)
for(ig in 1:ngrad){ 
z <- 0
for(ab in 1:Lind2) {
z <- z+4*pi*(-1)^(a[ab]/2)*Inn[a[ab]/2+1,n1,n2]*Qlmm[i1,i2,ab]*Yabm[ab,ig]
}
Kofg[i1,i2,n1,n2,ig] <- z
# use symmetry
Kofg[i2,i1,n1,n2,ig] <- z
Kofg[i2,i1,n2,n1,ig] <- z
Kofg[i1,i2,n2,n1,ig] <- z
}
Kofg
}
