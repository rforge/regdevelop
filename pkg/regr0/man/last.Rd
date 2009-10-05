\name{last}
\alias{last}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Last Elements of a Vector }
\description{
  Selects or drops the last element or the last \code{n} elements of a vector 
}
\usage{
last(data, n = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{vector from which to select or drop}
  \item{n}{if >0, \code{last} selects the last \code{n} elements form
    the result.\cr 
    if <0, the last \code{abs(n)} elements are dropped, and the first
    \code{length(data)-abs(n)} elements from the result}
}
\value{
  The selected elements of the vector
}
\author{Werner Stahel}
\note{
  This is a very simple function. It is defined mainly for selecting
  from the results of other functions without storing them.  
}
\examples{
  x <- runif(rpois(1,10))
  last(sort(x), 3)
  last(sort(x), -5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
