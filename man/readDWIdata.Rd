\name{readDWIdata}
\alias{readDWIdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readDWIdata(dirlist, format, nslice, gradient, order = NULL, xind = NULL, yind = NULL, zind = NULL, level = 0, mins0value = 0, maxvalue = 10000, voxelext = NULL, orientation = c(1, 3, 5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dirlist}{ ~~Describe \code{dirlist} here~~ }
  \item{format}{ ~~Describe \code{format} here~~ }
  \item{nslice}{ ~~Describe \code{nslice} here~~ }
  \item{gradient}{ ~~Describe \code{gradient} here~~ }
  \item{order}{ ~~Describe \code{order} here~~ }
  \item{xind}{ ~~Describe \code{xind} here~~ }
  \item{yind}{ ~~Describe \code{yind} here~~ }
  \item{zind}{ ~~Describe \code{zind} here~~ }
  \item{level}{ ~~Describe \code{level} here~~ }
  \item{mins0value}{ ~~Describe \code{mins0value} here~~ }
  \item{maxvalue}{ ~~Describe \code{maxvalue} here~~ }
  \item{voxelext}{ ~~Describe \code{voxelext} here~~ }
  \item{orientation}{ ~~Describe \code{orientation} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dirlist, format, nslice, gradient, order = NULL, xind = NULL, 
    yind = NULL, zind = NULL, level = 0, mins0value = 0, maxvalue = 10000, 
    voxelext = NULL, orientation = c(1, 3, 5)) 
{
    if (!(format \%in\% c("DICOM", "NIFTI", "ANALYZE", "AFNI"))) 
        stop("Cannot handle other formats then DICOM|NIFTI|ANALYZE|AFNI, found:", 
            format)
    if (any(sort((orientation)\%/\%2) != 0:2)) 
        stop("invalid orientation \n")
    if (dim(gradient)[2] == 3) 
        gradient <- t(gradient)
    if (dim(gradient)[1] != 3) 
        stop("Not a valid gradient matrix")
    ngrad <- dim(gradient)[2]
    s0ind <- (1:ngrad)[apply(abs(gradient), 2, max) == 0]
    if (is.null(zind)) 
        zind <- 1:nslice
    filelist <- NULL
    for (dd in dirlist) filelist <- c(filelist, paste(dd, list.files(dd), 
        sep = ""))
    if (format == "DICOM") {
        if (length(filelist) != ngrad * nslice) 
            stop("Number of found files does not match ngrad*nslice", 
                length(filelist))
        if (is.null(order)) {
            order <- 1:(ngrad * nslice)
        }
        else {
            if (length(order) != ngrad * nslice) 
                stop("Length of order vector does not match ngrad*nslice")
        }
        dim(order) <- c(nslice, ngrad)
        order <- order[zind, ]
        dim(order) <- NULL
        filelist <- filelist[order]
    }
    else {
        if (length(filelist) != ngrad) 
            stop("Number of found files does not match ngrad", 
                length(filelist), "\nPlease provide each gradient cube in a separate file.")
        if (is.null(order)) {
            order <- 1:ngrad
        }
        else {
            if (length(order) != ngrad) 
                stop("Length of order vector does not match ngrad")
        }
        filelist <- filelist[order]
    }
    cat("Start reading data", date(), proc.time(), "\n")
    si <- numeric()
    cat("\n")
    ddim <- NULL
    for (ff in filelist) {
        cat(".")
        if (format == "DICOM") {
            data <- read.DICOM(ff)
        }
        else if (format == "NIFTI") {
            data <- read.NIFTI(ff, setmask = FALSE)
        }
        else if (format == "ANALYZE") {
            data <- read.ANALYZE(ff, setmask = FALSE)
        }
        else if (format == "AFNI") {
            data <- read.AFNI(ff, setmask = FALSE)
        }
        if (is.null(ddim)) 
            ddim <- c(data$dim[1:2], nslice, ngrad)
        if (is.null(voxelext)) {
            if (!is.null(data$delta)) {
                if (!prod(voxelext == data$delta)) 
                  warning("Voxel extension", voxelext, "is not found in data:", 
                    data$delta)
                voxelext <- data$delta
            }
            else {
                warning("Voxel extension neither found nor given!")
            }
        }
        if (is.null(xind)) 
            xind <- 1:data$dim[1]
        if (is.null(yind)) 
            yind <- 1:data$dim[2]
        if (format == "DICOM") {
            si <- c(si, extract.data(data)[xind, yind])
        }
        else {
            ttt <- extract.data(data)
            si <- c(si, ttt[xind, yind, zind, ])
        }
    }
    cat("\n")
    dim(si) <- c(length(xind), length(yind), length(zind), ngrad)
    dimsi <- dim(si)
    cat("Data successfully read", date(), proc.time(), "\n")
    xyz <- (orientation)\%/\%2 + 1
    swap <- orientation - 2 * (orientation\%/\%2)
    if (any(xyz != 1:3)) {
        abc <- 1:3
        abc[xyz] <- abc
        si <- aperm(si, c(abc, 4))
        swap[xyz] <- swap
        voxelext[xyz] <- voxelext
        dimsi[xyz] <- dimsi[1:3]
        ddim[xyz] <- ddim[1:3]
        gradient[xyz, ] <- gradient
    }
    if (swap[1] == 1) {
        si <- si[dimsi[1]:1, , , ]
        gradient[1, ] <- -gradient[1, ]
    }
    if (swap[2] == 1) {
        si <- si[, dimsi[2]:1, , ]
        gradient[2, ] <- -gradient[2, ]
    }
    if (swap[3] == 0) {
        si <- si[, , dimsi[3]:1, ]
        gradient[3, ] <- -gradient[3, ]
    }
    si <- .Fortran("initdata", si = as.integer(si), as.integer(dimsi[1]), 
        as.integer(dimsi[2]), as.integer(dimsi[3]), as.integer(dimsi[4]), 
        as.integer(maxvalue), PACKAGE = "dti")$si
    dim(si) <- dimsi
    level <- max(mins0value, level * mean(si[, , , s0ind][si[, 
        , , s0ind] > 0]))
    ddim0 <- as.integer(ddim)
    ddim <- as.integer(dim(si)[1:3])
    cat("Create auxiliary statistics", date(), proc.time(), " \n")
    btb <- create.designmatrix.dti(gradient)
    rind <- replind(gradient)
    invisible(new("dtiData", si = si, btb = btb, ngrad = ngrad, 
        s0ind = s0ind, replind = rind, ddim = ddim, ddim0 = ddim0, 
        xind = xind, yind = yind, zind = zind, level = level, 
        voxelext = voxelext, orientation = as.integer(c(0, 2, 
            5)), source = paste(dirlist, collapse = "|")))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
