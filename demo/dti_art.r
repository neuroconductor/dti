lambda <- 47
sigma <- 1600
rho <- 1
ddim <- c(64,64,26)
ngrad <- 25
factor <- 2.5

bvec <- read.table(system.file("data/b-directions.txt",package="dti"))
bvec <- t(bvec)
btb <- matrix(0,6,ngrad+1)
btb[1,] <- bvec[1,]*bvec[1,]
btb[4,] <- bvec[2,]*bvec[2,]
btb[6,] <- bvec[3,]*bvec[3,]
btb[2,] <- 2*bvec[1,]*bvec[2,]
btb[3,] <- 2*bvec[1,]*bvec[3,]
btb[5,] <- 2*bvec[2,]*bvec[3,]

eta <- function(ai){
  aindex <- function(tensor){
    values <- eigen(matrix(tensor[c(1,2,3,2,4,5,3,5,6)],3,3))$values
    sqrt(3/2*sum((values-mean(values))^2)/sum(values^2))
  }
  risk <- function(par,ai) (ai-aindex((1-par)*c(1,0,0,0,0,0)+par*c(1,0,0,1,0,1)))^2
  optimize(f = risk, interval = c(0,1),ai=ai)$minimum
}

cat("create Phantom data\n")
etas <- numeric(1001)
for(i in 1:1001) etas[i] <- eta((i-1)/1000)

ind <- array(0,ddim)
dtiso <- array(c(1,0,0,1,0,1),dim=c(6,ddim))
phi <- (1:1000)*2*pi/1000 
rad1 <- 7
rad2 <- 8
for(rad in seq(rad1,rad2,.25)){
  x <- as.integer(rad*sin(phi)+32.5)
  y <- as.integer(rad*cos(phi)+32.5)
  for(i in 1:125) ind[x[i],y[i],2:25] <- 0
  for(i in 126:250) ind[x[i],y[i],2:25] <- .6
  for(i in 251:375) ind[x[i],y[i],2:25] <- .2
  for(i in 376:500) ind[x[i],y[i],2:25] <- .8
  for(i in 501:625) ind[x[i],y[i],2:25] <- .4
  for(i in 626:750) ind[x[i],y[i],2:25] <- .7
  for(i in 751:875) ind[x[i],y[i],2:25] <- .3
  for(i in 876:1000) ind[x[i],y[i],2:25] <- .9
  for(i in 1:1000){
    etai <- etas[ind[x[i],y[i],12]*1000+1]
    dtiso[,x[i],y[i],] <- (1-etai)*c(0,0,0,0,0,1)+etai*c(1,0,0,1,0,1)
  }
}

rad1 <- 13
rad2 <- 14
for(rad in seq(rad1,rad2,.25)){
  sphi <- sin(phi)
  cphi <- cos(phi)
  x <- as.integer(rad*sphi+32.5)
  y <- as.integer(rad*cphi+32.5)
  for(i in 1:1000) ind[x[i],y[i],2:4] <- .6
  for(i in 1:1000) ind[x[i],y[i],5:7] <- .3
  for(i in 1:1000) ind[x[i],y[i],8:10] <- .9
  for(i in 1:1000) ind[x[i],y[i],11:13] <- 0
  for(i in 1:1000) ind[x[i],y[i],14:16] <- .5
  for(i in 1:1000) ind[x[i],y[i],17:19] <- .2
  for(i in 1:1000) ind[x[i],y[i],20:22] <- .8
  for(i in 1:1000) ind[x[i],y[i],23:25] <- .4 
  for(j in 1:26) {
    etai <- etas[ind[x[i],y[i],j]*1000+1]
    for(i in 1:1000) dtiso[,x[i],y[i],j] <- (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],0,sphi[i]^2,0,0)+etai*c(1,0,0,1,0,1)
  }
}
rad1 <- 19
rad2 <- 20
for(rad in seq(rad1,rad2,.25)){
  sphi <- sin(phi)
  cphi <- cos(phi)
  x <- as.integer(rad*sphi+32.5)
  y <- as.integer(rad*cphi+32.5)
  for(i in 1:1000) ind[x[i],y[i],2:3] <- .6
  for(i in 1:1000) ind[x[i],y[i],4:5] <- .2
  for(i in 1:1000) ind[x[i],y[i],6:7] <- .8
  for(i in 1:1000) ind[x[i],y[i],8:9] <- 0
  for(i in 1:1000) ind[x[i],y[i],10:11] <- .5
  for(i in 1:1000) ind[x[i],y[i],12:13] <- .9
  for(i in 1:1000) ind[x[i],y[i],14:15] <- .2
  for(i in 1:1000) ind[x[i],y[i],16:17] <- .6 
  for(i in 1:1000) ind[x[i],y[i],18:19] <- 0 
  for(i in 1:1000) ind[x[i],y[i],20:21] <- .9 
  for(i in 1:1000) ind[x[i],y[i],22:23] <- .3 
  for(i in 1:1000) ind[x[i],y[i],24:25] <- .7
  for(j in 1:26) {
    etai <- etas[ind[x[i],y[i],j]*1000+1]
    for(i in 1:1000) dtiso[,x[i],y[i],j] <- (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],0,sphi[i]^2,0,0)+etai*c(1,0,0,1,0,1)
  }
}

#
#   to slow
#
rad1 <- 25
rad2 <- 26
for(rad in seq(rad1,rad2,.25)){
   x <- as.integer(rad*sin(phi)+32.5)
   y <- as.integer(rad*cos(phi)+32.5)
   sphi <- sin(phi)
   cphi <- cos(phi)
   for(i in 1:1000) ind[x[i],y[i],3:5] <- (1+sin(phi[i]))/3
   for(i in 1:1000) ind[x[i],y[i],6:8] <- (1+cos(phi[i]))/3
   for(i in 1:1000) ind[x[i],y[i],9:11] <- 0
   for(i in 1:1000) ind[x[i],y[i],12:13] <- (1+sin(phi[i]))/2
   for(i in 1:1000) ind[x[i],y[i],14:15] <- (1+cos(phi[i]))/2
   for(i in 1:1000) ind[x[i],y[i],16:18] <- 0
   for(i in 1:1000) ind[x[i],y[i],19:21] <- (2+sin(phi[i]))/3
   for(i in 1:1000) ind[x[i],y[i],22:24] <- (2+cos(phi[i]))/3
   for(j in 1:26) {
      for(i in 1:1000){
         etai <- etas[ind[x[i],y[i],j]*1000+1]
         dtiso[,x[i],y[i],j] <- (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],0,sphi[i]^2,0,0)+etai*c(1,0,0,1,0,1)
      }
     cat(".")
   }
}

for( i in 1:64) for (j in 1:64){
  if(max(ind[i,j,])==0&&((i-32.5)^2+(j-32.5)^2>26^2)) {
    dtiso[,i,j,]<-0
  }
}

dtiso <- factor*dtiso
#
#  yields a maximum eigenvalue of 2.5 within a cylinder of radius 26, zero tensor outside this cylinder
#
cat("show images slice 14 ")
par(mfrow=c(2,3),mar=c(1,1,1,.1))
for(i in 1:6) image(dtiso[i,,,14],col=grey((0:255)/255))

project.cylinder <- function(obj,radius,phi=(1:1000)*2*pi/1000){
  ddim <- dim(obj)
  img <- matrix(0,length(phi),ddim[3])
  x <- as.integer(radius*sin(phi)+32.5)
  y <- as.integer(radius*cos(phi)+32.5)
  for ( i in 1:length(phi)) img[i,] <- obj[x[i],y[i],]
  img
}

par(mfrow=c(2,3),mar=c(1,1,1,1))
image(ind[,,22],col=grey((0:255)/255))
image(project.cylinder(ind,7.5),col=grey((0:255)/255))
image(project.cylinder(ind,13),col=grey((0:255)/255))
image(project.cylinder(ind,19),col=grey((0:255)/255))
image(project.cylinder(ind,25),col=grey((0:255)/255))

s0offa <- read.table(system.file("data/S0ofFA.txt",package="dti"))
s0 <- s0offa[as.integer(as.vector(ind)*500+1),2]
dim(s0) <- dim(ind)
for( i in 1:64) for (j in 1:64){
  if(max(ind[i,j,])==0&&((i-32.5)^2+(j-32.5)^2>26^2)) s0[i,j,]<-0
}

createdata.dti2 <- function(file,dtensor,btb,s0,sigma,level=250){
  ngrad <- dim(btb)[2]
  ddim <- dim(s0)
  dim(dtensor)<-c(6,prod(ddim))
  dtensor <- t(dtensor)
  si <- exp(-dtensor%*%btb)*as.vector(s0)
  dim(si)<-c(ddim,ngrad)
  for (i in 1:ddim[3]) {
    s0[,,i] <- abs(fft(fft(s0[,,i])+complex(real=rnorm(64*64,0,sigma),imaginary=rnorm(64*64,0,sigma)),inverse=TRUE))/ddim[1]/ddim[2]
  }
  for (j in 1:ngrad) {
    for (i in 1:ddim[3]) {
      si[,,i,j] <- abs(fft(fft(si[,,i,j])+complex(real=rnorm(64*64,0,sigma),imaginary=rnorm(64*64,0,sigma)),inverse=TRUE))/ddim[1]/ddim[2]
    }
  }
  con <- file(file,"wb")
  writeBin(as.integer(si),con,2)
  close(con)
}

expected.dti <- function(file,dtensor,btb,s0,sigma,level=250){
  ngrad <- dim(btb)[2]
  ddim <- dim(s0)
  dim(dtensor)<-c(6,prod(ddim))
  dtensor <- t(dtensor)
  si <- exp(-dtensor%*%btb)*as.vector(s0)
  dim(si)<-c(ddim,ngrad)
  s0 <- s0+sigma^2/8100/s0
  si <- si+sigma^2/8100/si
#
#  1/ Bias S_i ~ 8100 / sigma^2 * true S_i    
#
  con <- file(file,"wb")
  writeBin(as.integer(si),con,2)
  close(con)
}


#
#   create phantom - object
#


dt0 <- new("dtiTensor",D=dtiso, th0= s0,sigma=array(0,dim(dtiso)[-1]), scorr=array(0,dim=c(5,5,3)), bw = c(0,0,0), mask=array(TRUE,dim=dim(dtiso)[-1]), ddim=dim(dtiso)[-1], ddim0=dim(dtiso)[-1], method="unknown")
dt0aniso <- dtiIndices(dt0)

sigma <- 0.01
expected.dti("S_noise_all",dtiso,btb,s0,sigma)
dtobj <- dtiData(bvec,paste("S_noise_all",sep=""),ddim)
dthat0 <- dtiTensor(dtobj)
dthat0aniso <- dtiIndices(dthat0)

cat("create noisy data\n")
set.seed(1)
sigma <- .01
createdata.dti2("S_noise_all",dtiso,btb,s0,sigma)

#
#  Bias of FA can not be avoided since the Expected S_0 S_b  and therefore the Expected Tensor 
#  differ from the "True" quantities 
#  We may be better off if we compare the FA with the FA of the Expected tensor computet from E S_0 and E S_b
#
#

# Read noisy data 
dtobj <- dtiData(bvec,paste("S_noise_all",sep=""),ddim)

# estimate noisy tensor too
cat("calculating noisy tensor\n")
dthat1 <- dtiTensor(dtobj)
dthat1aniso <- dtiIndices(dthat1)
cat("Adaptive smoothing \n")

dthat4 <- dti.smooth(dtobj,hmax=4,graph=TRUE,lambda=lambda,minanindex=0,slice=15,rho=rho,lseq=NULL)
dthat4aniso <- dtiIndices(dthat4)
par(mfrow=c(1,3))
plot(dt0aniso,slice=15)
plot(dthat1aniso,slice=15)
plot(dthat4aniso,slice=15)
ind2 <- logical(length(ind))
dim(ind2) <- dim(ind)
for( i in 1:64) for (j in 1:64){
   ind2[i,j,] <- ((i-32.5)^2+(j-32.5)^2<=26^2)
}

cat("\n")
cat("Parameters:","lambda=",lambda,"sigma=",sigma,"rho=",rho,"\n")
cat("MSE for FA","Original",mean((dt0aniso$fa-dthat1aniso$fa)[ind2]^2),"Smooth",mean((dt0aniso$fa-dthat4aniso$fa)[ind2]^2),"\n")
cat("Bias for FA","Original",mean((dt0aniso$fa-dthat1aniso$fa)[ind2]),"Smooth",mean((dt0aniso$fa-dthat4aniso$fa)[ind2]),"\n")
eigenv0 <- matrix(dt0aniso$eigenv[,,,,3],prod(dim(dthat1aniso$eigenv)[1:3]),3)
eigenv1 <- matrix(dthat1aniso$eigenv[,,,,3],prod(dim(dthat1aniso$eigenv)[1:3]),3)
eigenv4 <- matrix(dthat4aniso$eigenv[,,,,3],prod(dim(dthat1aniso$eigenv)[1:3]),3)
cat("Risk for Andir","Original",mean(dt0aniso$fa[ind2]*(1-abs(apply(eigenv0*eigenv1,1,sum))[ind2])),"Smooth",
mean(dt0aniso$fa[ind2]*(1-abs(apply(eigenv0*eigenv4,1,sum))[ind2])),"\n")

project2.cylinder <- function(obj,radius,phi=(1:1000)*2*pi/1000){
  ddim <- dim(obj)
  img <- array(0,c(length(phi),ddim[3],3))
  x <- as.integer(radius*sin(phi)+32.5)
  y <- as.integer(radius*cos(phi)+32.5)
  for ( i in 1:length(phi)) img[i,,] <- obj[x[i],y[i],,]
  img
}

projandir.image <- function(dtobject,minanindex=0,radius=10){
  if(!("dtiIndices" %in% class(dtobject))) stop("Not an dti-object")
  adimpro <- require(adimpro)
  anindex <- project.cylinder(dtobject$fa,radius=radius,phi=(1:256)*2*pi/256)
  andirection <- project2.cylinder(dtobject$eigenv[,,,,3],radius=radius,phi=(1:256)*2*pi/256)
  dimg <- dim(anindex)
  anindex[anindex>1]<-0
  anindex[anindex<0]<-0
  andirection <- aperm(andirection,c(3,1,2))
  andirection[1,,] <- abs(andirection[1,,])
  andirection[2,,] <- abs(andirection[2,,])
  andirection[3,,] <- abs(andirection[3,,])
  andirection <- aperm(andirection,c(2,3,1))
  andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
  dim(andirection)<-c(dimg,3)
  andirection <- andirection[,rep(1:26,rep(4,26)),]
  andirection <- make.image(andirection)
  invisible(andirection)
}
andir2.image <- function(dtobject,slice=1,method=1,quant=0,minanindex=NULL,show=TRUE,xind=NULL,yind=NULL,...){
  if(!("dtiIndices" %in% class(dtobject))) stop("Not an dtiIndices-object")
  adimpro <- require(adimpro)
  anindex <- dtobject$fa
  dimg <- dim(anindex)[1:2]
  if(is.null(xind)) xind <- 1:dimg[1]
  if(is.null(yind)) yind <- 1:dimg[2]
  if(is.null(slice)) slice <- 1
  anindex <- anindex[xind,yind,slice]
  dimg <- dim(anindex)[1:2]
  andirection <- dtobject$eigenv[xind,yind,slice,,3]
  anindex[anindex>1]<-0
  anindex[anindex<0]<-0
  dim(andirection)<-c(prod(dimg),3)
  andirection <- t(andirection)
  if(is.null(minanindex)) minanindex <- quantile(anindex,quant)
  if(method==1) {
    andirection[1,] <- abs(andirection[1,])
    andirection[2,] <- abs(andirection[2,])
    andirection[3,] <- abs(andirection[3,])
  } else {
    ind<-andirection[1,]<0
    andirection[,ind] <- - andirection[,ind]
    andirection[2,] <- (1+andirection[2,])/2
    andirection[3,] <- (1+andirection[3,])/2
  }
  andirection <- t(andirection)
  andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minanindex)
  dim(andirection)<-c(dimg,3)
  if(adimpro) {
    andirection <- andirection[rep(1:64,rep(4,64)),rep(1:64,rep(4,64)),]
    andirection <- make.image(andirection)
    if(show) show.image(andirection,...)
  } else if(show) {
    dim(anindex) <- dimg
    image(anindex,...)
  }
  invisible(andirection)
} 

img0s <- andir2.image(dt0aniso,slice=15,minanindex=0,show=FALSE)
img1s <- andir2.image(dthat1aniso,slice=15,minanindex=0,show=FALSE)
img4s <- andir2.image(dthat4aniso,slice=15,minanindex=0,show=FALSE)
eps <- .5
img0z1 <- projandir.image(dt0aniso,minanindex=0,radius=7.+eps)
img1z1 <- projandir.image(dthat1aniso,minanindex=0,radius=7.+eps)
img4z1 <- projandir.image(dthat4aniso,minanindex=0,radius=7.+eps)

img0z2 <- projandir.image(dt0aniso,minanindex=0,radius=13.+eps)
img1z2 <- projandir.image(dthat1aniso,minanindex=0,radius=13.+eps)
img4z2 <- projandir.image(dthat4aniso,minanindex=0,radius=13.+eps)

img0z3 <- projandir.image(dt0aniso,minanindex=0,radius=19.+eps)
img1z3 <- projandir.image(dthat1aniso,minanindex=0,radius=19.+eps)
img4z3 <- projandir.image(dthat4aniso,minanindex=0,radius=19.+eps)

img0z4 <- projandir.image(dt0aniso,minanindex=0,radius=25.+eps)
img1z4 <- projandir.image(dthat1aniso,minanindex=0,radius=25.+eps)
img4z4 <- projandir.image(dthat4aniso,minanindex=0,radius=25.+eps)
write.image(img0s,"img0sa.png",depth=8)
write.image(img1s,"img1sa.png",depth=8)
write.image(img4s,"img4sa.png",depth=8)

write.image(img0z1,"img0z1a.png",depth=8)
write.image(img1z1,"img1z1a.png",depth=8)
write.image(img4z1,"img4z1a.png",depth=8)

write.image(img0z2,"img0z2a.png",depth=8)
write.image(img1z2,"img1z2a.png",depth=8)
write.image(img4z2,"img4z2a.png",depth=8)

write.image(img0z3,"img0z3a.png",depth=8)
write.image(img1z3,"img1z3a.png",depth=8)
write.image(img4z3,"img4z3a.png",depth=8)

write.image(img0z4,"img0z4a.png",depth=8)
write.image(img1z4,"img1z4a.png",depth=8)
write.image(img4z4,"img4z4a.png",depth=8)
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,248 \"a)\"' img0sa.png img0s.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,248 \"f)\"' img1sa.png img1s.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,248 \"k)\"' img4sa.png img4s.png")

system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"b)\"' img0z1a.png img0z1.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"g)\"' img1z1a.png img1z1.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"l)\"' img4z1a.png img4z1.png")

system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"c)\"' img0z2a.png img0z2.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"h)\"' img1z2a.png img1z2.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"m)\"' img4z2a.png img4z2.png")

system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"d)\"' img0z3a.png img0z3.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"i)\"' img1z3a.png img1z3.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"n)\"' img4z3a.png img4z3.png")

system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"e)\"' img0z4a.png img0z4.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"j)\"' img1z4a.png img1z4.png")
system("convert -fill white -font Courier-Bold -pointsize 24 -draw 'text 6,96 \"o)\"' img4z4a.png img4z4.png")

system("convert -border 4x0 img0s.png img1s.png img4s.png +append art_s.png")
system("convert -border 4x0  img0z1.png img1z1.png img4z1.png +append art_z1.png")
system("convert -border 4x0  img0z2.png img1z2.png img4z2.png +append art_z2.png")
system("convert -border 4x0  img0z3.png img1z3.png img4z3.png +append art_z3.png")
system("convert -border 4x0  img0z4.png img1z4.png img4z4.png +append art_z4.png")
system(paste("convert -border 0x4  art_s.png art_z1.png art_z2.png art_z3.png art_z4.png -append art_all_",as.integer(sigma),"_",as.integer(lambda),"_",signif(rho,2),".png",sep=""))

system(paste("cp art_all_",as.integer(sigma),"_",as.integer(lambda),"_",signif(rho,2),".png ../neuroimage-revision1/fig2.png",sep=""))
ind1 <- ind

for( i in 1:64) for (j in 1:64){
   if(((i-32.5)^2+(j-32.5)^2>562)) ind1[i,j,]<- -1
}

ind1[ind1==0.2] <- 1
ind1[ind1==0.3] <- 2
ind1[ind1==0.4] <- 3
ind1[ind1==0.5] <- 4
ind1[ind1==0.6] <- 5
ind1[ind1==0.7] <- 6
ind1[ind1==0.8] <- 7
ind1[ind1==0.9] <- 8

ae01 <- sb01 <- ae1 <- sb1 <- list(NULL)
for(i in 0:8) ae1[[i+1]]<-as.vector(abs(dt0aniso$fa-dthat1aniso$fa)[ind1==i])
for(i in 0:8) sb1[[i+1]]<-as.vector(-(dt0aniso$fa-dthat1aniso$fa)[ind1==i])
for(i in 0:8) ae01[[i+1]]<-as.vector(abs(dthat0aniso$fa-dthat1aniso$fa)[ind1==i])
for(i in 0:8) sb01[[i+1]]<-as.vector(-(dthat0aniso$fa-dthat1aniso$fa)[ind1==i])

ae04 <- sb04 <- ae4 <- sb4 <- list(NULL)
for(i in 0:8) ae4[[i+1]]<-as.vector(abs(dt0aniso$fa-dthat4aniso$fa)[ind1==i])
for(i in 0:8) sb4[[i+1]]<-as.vector(-(dt0aniso$fa-dthat4aniso$fa)[ind1==i])
for(i in 0:8) ae04[[i+1]]<-as.vector(abs(dthat0aniso$fa-dthat4aniso$fa)[ind1==i])
for(i in 0:8) sb04[[i+1]]<-as.vector(-(dthat0aniso$fa-dthat4aniso$fa)[ind1==i])

names(ae1) <- names(sb1) <- names(ae4) <- names(sb4) <-
names(ae01) <- names(sb01) <- names(ae04) <- names(sb04) <- 
c("FA=0","FA=0.2","FA=0.3","FA=0.4",
"FA=0.5","FA=0.6","FA=0.7","FA=0.8","FA=0.9")

png("anindex.png",width=1000,height=800)
par(mfrow=c(2,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(ae1,notch=TRUE,outline=FALSE,ylim=c(0,.3),xlab="FA - value",ylab="Absolute error FA")
title("Voxelwise estimates")
boxplot(ae4,notch=TRUE,outline=FALSE,ylim=c(0,.3),xlab="FA - value",ylab="Absolute error FA")
title("Smoothed estimates")
boxplot(sb1,notch=TRUE,outline=FALSE,ylim=c(-.32,.25),xlab="FA - value",ylab="Bias FA")
title("Voxelwise estimates")
boxplot(sb4,notch=TRUE,outline=FALSE,ylim=c(-.32,.25),xlab="FA - value",ylab="Bias FA")
title("Smoothed  estimates")
dev.off()
system("cp anindex.png ../fig4.png")
pdf("anindex.pdf",width=1000/72,height=800/72)
par(mfrow=c(2,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(ae1,notch=TRUE,outline=FALSE,ylim=c(0,.3),xlab="FA - value",ylab="Absolute error FA")
title("Voxelwise estimates")
boxplot(ae4,notch=TRUE,outline=FALSE,ylim=c(0,.3),xlab="FA - value",ylab="Absolute error FA ")
title("Smoothed estimates")
boxplot(sb1,notch=TRUE,outline=FALSE,ylim=c(-.32,.25),xlab="FA - value",ylab="Error of FA estimate")
title("Voxelwise estimates")
boxplot(sb4,notch=TRUE,outline=FALSE,ylim=c(-.32,.25),xlab="FA - value",ylab="Error of FA estimate")
title("Smoothed  estimates")
dev.off()
system("cp anindex.pdf ../fig4.pdf")
system("cp anindex.png ../fig4.png")

pdf("anindexb.pdf",width=1000/72,height=800/72)
par(mfrow=c(2,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(ae01,notch=TRUE,outline=FALSE,ylim=c(0,.26),xlab="FA - value",ylab="Absolute error FA")
title("(a) Voxelwise estimates")
boxplot(ae04,notch=TRUE,outline=FALSE,ylim=c(0,.26),xlab="FA - value",ylab="Absolute error FA")
title("(b) Smoothed estimates")
boxplot(sb01,notch=TRUE,outline=FALSE,ylim=c(-.24,.3),xlab="FA - value",ylab="Error of FA estimate")
title("(c) Voxelwise estimates")
boxplot(sb04,notch=TRUE,outline=FALSE,ylim=c(-.24,.3),xlab="FA - value",ylab="Error of FA estimate")
title("(d) Smoothed  estimates")
dev.off()
system("cp anindexb.pdf ../fig4b.pdf")

pdf("anindexc.pdf",width=1050/72,height=470/72)
par(mfrow=c(1,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(sb01,notch=TRUE,outline=FALSE,ylim=c(-.24,.3),xlab="FA - value",ylab="Error of FA estimate")
title("(a) Voxelwise estimates")
boxplot(sb04,notch=TRUE,outline=FALSE,ylim=c(-.24,.3),xlab="FA - value",ylab="Error of FA estimate")
title("(b) Smoothed  estimates")
dev.off()
system("cp anindexc.pdf ../fig4c.pdf")

ad1 <- ad4 <- ad01 <- ad04 <- list(NULL)
eigenv0<-matrix(dt0aniso$eigenv[,,,,3],prod(dim(dthat1aniso$eigenv)[1:3]),3)
eigenv0b<-matrix(dthat0aniso$eigenv[,,,,3],prod(dim(dthat1aniso$eigenv)[1:3]),3)
eigenv1<-matrix(dthat1aniso$eigenv[,,,,3],prod(dim(dthat1aniso$eigenv)[1:3]),3)
eigenv4<-matrix(dthat4aniso$eigenv[,,,,3],prod(dim(dthat1aniso$eigenv)[1:3]),3)
ind1 <- as.vector(ind1)
edir1 <- acos(abs(apply(eigenv0*eigenv1,1,sum))*(1-1e-12))
edir4 <- acos(abs(apply(eigenv0*eigenv4,1,sum))*(1-1e-12))
edir01 <- acos(abs(apply(eigenv0b*eigenv1,1,sum))*(1-1e-12))
edir04 <- acos(abs(apply(eigenv0b*eigenv4,1,sum))*(1-1e-12))
#edir1 <- 1-abs(apply(eigenv0*eigenv1,1,sum))
#edir4 <- 1-abs(apply(eigenv0*eigenv4,1,sum))
for(i in 1:8) ad1[[i]]<-edir1[ind1==i]
for(i in 1:8) ad4[[i]]<-edir4[ind1==i]
for(i in 1:8) ad01[[i]]<-edir01[ind1==i]
for(i in 1:8) ad04[[i]]<-edir04[ind1==i]
names(ad1) <- names(ad4) <-names(ad01) <- names(ad04) <-
 c("FA=0.2","FA=0.3","FA=0.4","FA=0.5","FA=0.6","FA=0.7","FA=0.8","FA=0.9")

png("andir.png",width=1000,height=450)
par(mfrow=c(1,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(ad1,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("Voxelwise estimates")
boxplot(ad4,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("Smoothed estimates")
dev.off()

system("cp andir.png ../fig3.png")
pdf("andir.pdf",width=1000/72,height=450/72)
par(mfrow=c(1,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(ad1,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("(a) Voxelwise estimates")
boxplot(ad4,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("(b) Smoothed estimates")
dev.off()
system("cp andir.pdf ../fig3.pdf")

pdf("andirb.pdf",width=1050/72,height=470/72)
par(mfrow=c(1,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(ad01,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("(a) Voxelwise estimates")
boxplot(ad04,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("(b) Smoothed estimates")
dev.off()
system("cp andirb.pdf ../fig3b.pdf")

cat("MAE reduction",(unlist(lapply(ae01,mean))-unlist(lapply(ae04,mean)))/unlist(lapply(ae01,mean)),"\n")
cat("Bias reduction",(unlist(lapply(sb01,mean))-unlist(lapply(sb04,mean)))/unlist(lapply(sb01,mean)),"\n")
cat("radiant error reduction",(unlist(lapply(ad01,mean))-unlist(lapply(ad04,mean)))/unlist(lapply(ad01,mean)),"\n")



pdf("andir.pdf",width=1050/72,height=550/72)
par(mfrow=c(1,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(ad01,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("Voxelwise estimates")
boxplot(ad04,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("Smoothed estimates")
dev.off()
pdf("anindex.pdf",width=1050/72,height=550/72)
par(mfrow=c(1,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(sb01,notch=TRUE,outline=FALSE,ylim=c(-.24,.3),xlab="FA - value",ylab="Error of FA estimate")
title("Voxelwise estimates")
boxplot(sb04,notch=TRUE,outline=FALSE,ylim=c(-.24,.3),xlab="FA - value",ylab="Error of FA estimate")
title("Smoothed  estimates")
dev.off()
pdf("andir2.pdf",width=1050/72,height=525/72)
par(mfrow=c(1,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(ad01,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("Voxelwise estimates")
boxplot(ad04,notch=TRUE,outline=TRUE,ylim=c(0,0.7),xlab="FA - value",ylab="Absolute error direction (radian)")
title("Smoothed estimates")
dev.off()
pdf("anindex2.pdf",width=1050/72,height=525/72)
par(mfrow=c(1,2),mar=c(3.5,3.5,3,.25),mgp=c(2,1,0))
boxplot(sb01,notch=TRUE,outline=FALSE,ylim=c(-.24,.3),xlab="FA - value",ylab="Error of FA estimate")
title("Voxelwise estimates")
boxplot(sb04,notch=TRUE,outline=FALSE,ylim=c(-.24,.3),xlab="FA - value",ylab="Error of FA estimate")
title("Smoothed  estimates")
dev.off()
