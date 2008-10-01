\name{extract}
\alias{extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract information from objects of class dtiData, dtiTensor or dtiIndices }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract(obj, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ ~~Describe \code{obj} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
structure(function (obj, ...) 
standardGeneric("extract"), generic = structure("extract", package = "dti"), package = "dti", group = list(), valueClass = character(0), signature = "obj", default = <S4 object of class structure("MethodsList", package = "methods")>, skeleton = structure(function (obj, 
    ...) 
cat("Data extraction not defined for this class:", class(object), 
    "\n"), target = structure(character(0), .Names = character(0), class = structure("signature", package = "methods")), defined = structure(character(0), .Names = character(0), class = structure("signature", package = "methods")), generic = structure("extract", package = "dti"), class = structure("derivedDefaultMethod", package = "methods"))(obj, 
    ...), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
