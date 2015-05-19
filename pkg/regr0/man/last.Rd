\name{last}
\alias{last}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Last Elements of a Vector or of a Matrix}
\description{
  Selects or drops the last element or the last \code{n} elements of a
  vector or the last \code{n} rows or \code{ncol} columns of a matrix 
}
\usage{
last(data, n = NULL, ncol=NULL, drop=is.matrix(data))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{vector or matrix or data.frame from which to select or drop}
  \item{n}{if >0, \code{last} selects the last \code{n} elements (rows) form
    the result.\cr 
    if <0, the last \code{abs(n)} elements (rows) are dropped, and the first
    \code{length(data)-abs(n)} ones from the result}
  \item{ncol}{if \code{data} is a matrix or data.frame, the last
    \code{ncol} columns are selected (if \code{ncol} is positive)
    or dropped (if negative).}
  \item{drop}{if only one row or column of a matrix (or one column of a
    data.frame) is selected or left over, should the result be a vector
    or a row or column matrix (or one variable data.frame)}
}
\value{
  The selected elements of the vector or matrix or data.frame
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
##
  df <- data.frame(X=c(2,5,3,8), F=LETTERS[1:4], G=c(TRUE,FALSE,FALSE,TRUE))
  last(df,3,-2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
