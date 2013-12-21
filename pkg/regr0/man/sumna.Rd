\name{sumna}
\alias{sumna}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Count NAs}
\description{
Count the missing or non-finite values for each column of a matrix or
data.frame 
}
\usage{
sumna(object, inf = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a vector, matrix, or data.frame}
  \item{inf}{if TRUE, Inf and NaN values are counted along with NAs}
}
\value{
numerical vector containing the missing value counts for each column
}
\author{Werner A. Stahel, ETH Zurich}
\note{This is a simple shortcut for \code{apply(is.na(object),2,sum)}
  or \code{apply(!is.finite(object),2,sum)}
}
\seealso{\code{\link{is.na}}, \code{\link{is.finite}}}
\examples{
t.d <- data.frame(V1=c(1,2,NA,4), V2=c(11,12,13,Inf), V3=c(21,NA,23,Inf))
sumna(t.d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{NA}
