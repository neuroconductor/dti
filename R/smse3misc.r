expm <- function(m,eps=1e-12){
n <- dim(m)[1]
mpot <- diag(n)
ex <- diag(n)
jfac <- 1
j <- 1
while(max(abs(mpot/jfac/j))>eps&&j<36){
jfac <- jfac*j
mpot <- mpot%*%m
ex <- ex +mpot/jfac
j <- j+1
}
ex
}
expm <- function(m){
matrix(.Fortran("expm3", as.double(m),ex=double(9),
                    DUPL=FALSE,
                    PACKAGE="dti")$ex,3,3)
}
abofg <- function(g,eps=1e-7){
g <- g/sqrt(sum(g^2))
beta <- asin(g[1])
if(abs(g[1])>eps){
z <- g[3]/cos(beta)
if(abs(z) < 1) { 
gamma <- asin(z)
} else {
gamma <- sign(z)* pi/2
}
} else {
gamma <- 0
# any value is correct
}
c(beta,gamma)
}
abofg <- function(g){
dg <- dim(g)
ngrad <- if(!is.null(dg)) dg[2] else 1
bg <- .Fortran("abofg",as.double(g),
                    as.integer(ngrad),
                    bg=double(2*ngrad),
                    DUPL=FALSE,
                    PACKAGE="dti")$bg
dim(bg) <- c(2,ngrad)
bg
}

betagamma <- function(g1,g2){
bg1 <- abofg(g1)
bg2 <- abofg(g2)
b1 <- bg1[1]
b2 <- bg2[1]
dg <- bg1[2]-bg2[2]
z1 <- cos(b2)*sin(b1)-sin(b2)*cos(b1)*cos(dg)
bhat <- asin(z1)
z2 <- cos(b1)*sin(dg)/cos(bhat)
z2 <- sign(z2)*pmin(1,abs(z2))
ghat <- asin(z2)
if(any(is.nan(ghat))){
cat("b2",b2,"dg",dg,"\n")
print(rbind(b1,bhat,cos(b1)*sin(dg),cos(bhat),z2)[,is.na(asin(z2))])
ghat[is.nan(ghat)] <- pi/2
}
c(bhat,ghat)
}
betagamma <- function(g){
dg <- dim(g)
ngrad <- if(!is.null(dg)) dg[2] else 1
z <- .Fortran("bgstats",as.double(g),
                   as.integer(ngrad),
                   bg=double(2*ngrad),
                   bghat=double(2*ngrad*ngrad),
                   nbg=double(9*ngrad),
                   nbghat=double(9*ngrad*ngrad),
                   DUPL=FALSE,
                   PACKAGE="dti")[c("bg","bghat","nbg","nbghat")]
dim(z$bg) <- c(2,ngrad)
# sphaerische Coordinaten fuer Gradienten
dim(z$bghat) <- c(2,ngrad,ngrad)
# sphaerische Coordinaten fuer Gradienten-Paare
dim(z$nbg) <- c(3,3,ngrad)
# normalen-vektoren n1,n2,n3 Gradienten
dim(z$nbghat) <- c(3,3,ngrad,ngrad)
# normalen-vektoren n1,n2,n3 Gradienten-Paare
z
}
matrm <- function(b,g){
matrix(c(cos(b), 0, sin(b), sin(b)*sin(g), 
  cos(g), -cos(b) *sin(g), -cos(g)*sin(b), sin(g), cos(b)*cos(g)),3,3)
}
matrn4 <- function(b){
matrix(c(0,tan(b),0,-tan(b),0,-1/cos(b),0,1/cos(b),0),3,3)
}
getkappas <- function(grad,trace=0){
krit <- function(par,matm,m4,m5,m6){
#sum((matm-expm(par[1]*m4)%*%expm(par[2]*m5)%*%expm(par[3]*m6))^2)
.Fortran("k456krit",as.double(par),
                    as.double(matm),
                    as.double(m4),
                    as.double(m5),
                    as.double(m6),
                    erg=double(1),
                    DUPL=FALSE,
                    PACKAGE="dti")$erg
}
prta <- Sys.time()
cat("Start computing spherical distances",format(Sys.time()),"\n")
ngrad <- dim(grad)[2]
m5 <- matrix(c(0,0,1,0,0,0,-1,0,0),3,3)
m6 <- matrix(c(0,-1,0,1,0,0,0,0,0),3,3)
kappa456 <- array(0,c(3,ngrad,ngrad))
zbg <- betagamma(grad)
for(i in 1:ngrad) for(j in 1:ngrad) {
bg <- zbg$bghat[,i,j]
m4 <- matrn4(bg[1])
matm <- matrm(bg[1],bg[2])
k456 <- runif(3,-.01,.01)
z <-  optim(k456,krit,method="BFGS",matm=matm,m4=m4,m5=m5,m6=m6,control=list(trace=trace,reltol=1e-12,abstol=1e-12))
while(z$value>1e-8) {
#cat("i",i,"j",j,"value",z$value,"par",z$par,"\n")
k456 <- runif(3,-.01,.01)
z <- optim(k456,krit,method="BFGS",matm=matm,m4=m4,m5=m5,m6=m6,control=list(trace=trace,reltol=1e-12,abstol=1e-12))
#cat(" new value",z$value,"par",z$par,"\n")
}
kappa456[,i,j] <- z$par
}
while(any(abs(kappa456[2:3,,])>pi)){
kappa456[2:3,,][kappa456[2:3,,]< -pi] <- kappa456[2:3,,][kappa456[2:3,,]< -pi]+2*pi
kappa456[2:3,,][kappa456[2:3,,]> pi] <- kappa456[2:3,,][kappa456[2:3,,]> pi]-2*pi
}
prtb <- Sys.time()
cat("End computing spherical distances",format(Sys.time()),"\n")
list(k456=kappa456,bg=zbg$bg,bghat=zbg$bghat,nbg=zbg$nbg,nbghat=zbg$nbghat)
}

mthomogen <- function(object,minw=.1,maxangle=30){
andir <- extract(object,"andir")$andir
mix <- extract(object,"mix")$mix
order <- extract(object,"order")$order
mask <- extract(object,"mask")$mask
ddim <- object@ddim
z <- .Fortran("mthomog",
              as.double(andir),
              mix=as.double(mix),
              order=as.integer(order),
              as.integer(ddim[1]),
              as.integer(ddim[2]),
              as.integer(ddim[3]),
              as.integer(dim(mix)[1]),
              as.logical(mask),
              as.double(minw),
              as.double(maxangle/180*pi),
              as.double(object@voxelext),
              andir=as.double(andir),
              DUPL=TRUE,
              PACKAGE="dti")[c("andir","order","mix")]
object@orient <- array(.Fortran("parofor",
                                as.double(z$andir),
                                as.integer(prod(dim(mix))),
                                orient=double(2*prod(dim(mix))),
                                DUPL=FALSE,
                                PACKAGE="dti")$orient,c(2,dim(mix)))
object@mix <- array(z$mix,dim(mix))
object@order <- array(z$order,ddim)
object
}




