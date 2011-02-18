getkappa <- function(grad,dist){
kappa <- 1.25
ngrad <- dim(grad)[2]
w<-lkernse3(.9999,kappa,grad,c(1,1),dist)$w
vr <- sum(w)^2/sum(w^2)/ngrad
cat("kappa",kappa,"vr",vr,"\n")
while(vr>1+1e-4){
if(dist=="SE3") kappa <- kappa/1.25 else kappa <- kappa*1.25
w<-lkernse3(.9999,kappa,grad,c(1,1),dist)$w
vr <- sum(w)^2/sum(w^2)/ngrad
cat("kappa",kappa,"vr",vr,"\n")
}
while(vr<1+1e-4){
if(dist=="SE3") kappa <- kappa*1.01 else kappa <- kappa/1.01
w<-lkernse3(.9999,kappa,grad,c(1,1),dist)$w
vr <- sum(w)^2/sum(w^2)/ngrad
cat("kappa",kappa,"vr",vr,"\n")
}
kappa
}

gethseqse3 <- function(kstar,grad,vext=c(1,1),dist="SE3"){
h <- numeric(kstar)
ngrad <- dim(grad)[2]
cat("ngrad:",ngrad,"\n")
if(dist=="SE3"){
kappa <- 1/sqrt(.25+.17*ngrad)
} else {
kappa <- .28+.17*ngrad
}
#kappa <- getkappa(grad)
# thats korrect for optimized gradients
hakt <- 1.0
prt0 <- Sys.time()
for(k in 1:kstar){
w<-lkernse3(hakt,kappa,grad,vext,dist=dist)$w
vred <- sum(w)^2/sum(w^2)/ngrad/1.25^k
while(vred<1){
hakt <- hakt*1.01
w<-lkernse3(hakt,kappa,grad,vext,dist=dist)$w
vred <- sum(w)^2/sum(w^2)/ngrad/1.25^k
}
h[k]<-hakt
cat("k:",k,"vred:",signif(1.25^k,3),"h_k:",h[k],"sum(w)",
  signif(sum(w)/ngrad,3),"#w",signif(length(w)/ngrad,3),
  "time elapsed:",format(difftime(Sys.time(),prt0),digits=3),"\n")
}
h
}

