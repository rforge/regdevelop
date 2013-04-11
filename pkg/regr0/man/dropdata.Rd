\name{dropdata}
\alias{dropdata}
%\alias{seldata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Drop Observations from a Data.frame}
\description{
  Allows for dropping observations (rows) determined by row names or
  factor levels from a data.frame or matrix.  
}
\usage{
dropdata(data, rowid = NULL, incol = "row.names", colid = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data.frame of matrix}
  \item{rowid}{vector of character strings identifying the rows to be
    dropped}
  \item{incol}{name or index of the column used to identify the
    observations (rows)}
  \item{colid}{vector of character strings identifying the columns to be
    dropped}
}
\value{
  The data.frame or matrix without the dropped observations and/or
  variables. 
  Attributes are passed on.
}
\author{Werner A. Stahel, ETH Zurich}
\note{Ordinary subsetting by \code{[...,...]} drops attributes like
  \code{\link{doc}} or \code{\link{tit}}. Furthermore, the convenient
  way to drop rows or columns by giving negative indices to
  \code{[...,...]} cannot be used with names of rows or columns.
}
\seealso{
  \code{\link{subset}}
}
\examples{
dd <- data.frame(rbind(a=1:3,b=4:6,c=7:9,d=10:12))
dropdata(dd,"b")
dropdata(dd, col="X3")

d1 <- dropdata(dd,"d")
d2 <- dropdata(d1,"b")
naresid(attr(d2,"na.action"),as.matrix(d2))

dropdata(letters, 3:5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
