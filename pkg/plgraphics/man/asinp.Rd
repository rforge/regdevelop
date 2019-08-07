\name{asinp}
\alias{asinp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{arc sine Transformation}
\description{
  Calculates the sqrt arc sine of x/100, rescaled to be in the unit
  interval.\cr
  This transformation is useful for analyzing percentages or proportions
  of any kind.
}
\usage{
asinp(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of data values}
}
\value{
vector of transformed values}
\author{Werner A. Stahel, ETH Zurich}
\note{This very simple function is provided in order to simplify
  formulas. It has an attribute \code{"inverse"} that contains
  the inverse function, see example.  
}
\examples{
asinp(seq(0,100,10))
( y <- asinp(c(1,50,90,95,99)) )
attr(asinp, "inverse")(y)
}
\keyword{math}
