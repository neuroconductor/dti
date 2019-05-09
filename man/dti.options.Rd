\name{dti.options}
\alias{dti.options}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Set and manipulate image orientations for plots.
}
\description{
   The function can be used to adjust to radiological conventions
    in image displays.
}
\usage{
dti.options(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
  The following parameters can be used to determine the
  behaviour of the plot method for 3D image data in subsequent
  calls:
  \itemize{
     \item swapx - swap image x axis for display
     \item swapy - swap image y axis for display
     \item swapz - swap image z axis for display
}
all default to FALSE.
}
}
\value{
returns specified display orientations.
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }% use one of  RShowDoc("KEYWORDS")
