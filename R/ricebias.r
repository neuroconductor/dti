dwiRiceBias <- function(object,  ...) cat("No Rice Bias correction defined for this class:",class(object),"\n")

setGeneric("dwiRiceBias", function(object,  ...) standardGeneric("dwiRiceBias"))

setMethod("dwiRiceBias","dtiData",function(object, 
sigma=NULL, method="1stMoment") {
args <- object@call
if(is.null(sigma)||sigma<1){
cat("please provide a value for sigma \n Returning the original object\n")
return(object)
}
corrected <- FALSE
for(i in 1:length(args)){
if(length(grep("dwiRiceBias",args[i][[1]]))>0){
corrected <- TRUE
cat("Rice bais correction already performed by\n")
print(args[i][[1]])
cat("\n Returning the original object\n")
}
}
if(!corrected){
object@si <- array(ricebiascorr(object@si,sigma),dim(object@si))
object@call <- c(args,sys.call(-1))
}
invisible(object)
})



rrician <- function(x,m=0,s=1){
x1 <- rnorm(x,m,s)
x2 <- rnorm(x,0,s)
sqrt(x1*x1+x2*x2)
}

ricebiascorr <- function(x,s=1,method=1){
Lhalf <- function(x){
(1-x)*besselI(-x/2,0,TRUE) - x*besselI(-x/2,1,TRUE)
}
fofx <- function(x) sqrt(pi/2)*Lhalf(-x^2/2)
xt <- x/s
y <- xt
if(method==1){
xta <- 0
ind <- abs(xt-xta)>1e-6
while(sum(ind)>0){
xta <- xt
xt[ind] <- pmax(0,xt[ind]+y[ind]-fofx(xt[ind]))
ind <- abs(xt-xta)>1e-6
}
} else {
xt <- sqrt(pmax(0,y^2-2))
}
xt*s
}

