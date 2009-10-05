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
dropdata(data, rowid, incol = "row.names")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data.frame of matrix}
  \item{rowid}{vector of character strings identifying the rows to be
    dropped}
  \item{incol}{name or index of the column used to identify the
    observations (rows)}
}
\value{
  The data.frame or matrix without the dropped observations.
  Attributes are kept.
}
\author{Werner A. Stahel, ETH Zurich}
\note{Ordinary subsetting by \code{[...,...]} drops attributes like
  \code{\link{doc}} or \code{\link{tit}}. Furthermore, row or column
  names cannot be used like indices, simply by a preceding minus sign to
  drop the respective data.
}
\examples{
dd <- rbind(a=1:3,b=4:6,c=7:9)
dropdata(dd,"b")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
