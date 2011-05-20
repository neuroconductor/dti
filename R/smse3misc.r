expm <- function(m,eps=1e-12){
n <- dim(m)[1]
mpot <- diag(n)
ex <- diag(n)
jfac <- 1
j <- 1
while(max(abs(mpot/jfac/j))>eps){
jfac <- jfac*j
mpot <- mpot%*%m
ex <- ex +mpot/jfac
j <- j+1
}
ex
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
betagamma <- function(g1,g2){
bg1 <- abofg(g1)
bg2 <- abofg(g2)
cat("bg1",bg1,"bg2",bg2,"\n")
b1 <- bg1[1]
b2 <- bg2[1]
dg <- bg1[2]-bg2[2]
z1 <- cos(b2)*sin(b1)-sin(b2)*cos(b1)*cos(dg)
bhat <- asin(z1)
z2 <- cos(b1)*sin(dg)/cos(bhat)
if(any(is.nan(asin(z2)))){
cat("b2",b2,"dg",dg,"\n")
print(rbind(b1,bhat,cos(b1)*sin(dg),cos(bhat),z2)[,is.na(asin(z2))])
z2[is.nan(z2)] <- pi/2
}
ghat <- asin(z2)
c(bhat,ghat)
}
matrm <- function(b,g){
matrix(c(cos(b), 0, sin(b), sin(b)*sin(g), 
  cos(g), -cos(b) *sin(g), -cos(g)*sin(b), sin(g), cos(b)*cos(g)),3,3)
}
matrn4 <- function(b,g){
matrix(c(0,tan(b),0,-tan(b),0,-1/cos(b),0,1/cos(b),0),3,3)
}
getk456 <- function(g1,g2,trace=0){
m5 <- matrix(c(0,0,1,0,0,0,-1,0,0),3,3)
m6 <- matrix(c(0,-1,0,1,0,0,0,0,0),3,3)
bg <- betagamma(g1,g2)
m4 <- matrn4(bg[1],bg[2])
matm <- matrm(bg[1],bg[2]) 
cat("bg",bg,"\n")
k456 <- rep(0,3)
krit <- function(par,matm,m4,m5,m6){
sum((matm-expm(par[1]*m4)%*%expm(par[2]*m5)%*%expm(par[3]*m6))^2)
}
z <- optim(k456,krit,method="BFGS",matm=matm,m4=m4,m5=m5,m6=m6,control=list(trace=trace,reltol=1e-12,abstol=1e-12))
z
}



