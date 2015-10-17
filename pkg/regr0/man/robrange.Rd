\name{robrange}
\alias{robrange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust Range of Data }
\description{
  Determines a robust range of the data on the basis of the trimmed mean
  and mean absolute deviation
}
\usage{
robrange(data, trim = 0.2, fac = 3, na.rm=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ a vector of data. Missing values are dropped }
  \item{trim}{ trimming proportion }
  \item{fac}{ factor used for expanding the range, see Details }
  \item{na.rm}{logical: should NAs be removed? If FALSE, result will be
    NA if there are NAs in 'data'.}
}
\details{
  The function determines the trimmed mean \code{m} and then the "upper
  trimmed mean" \code{s} of absolute deviations from m, multiplied by
  \code{fac}. The robust minimum is then defined as \code{m-fac*s} or
  \code{min(data)}, whichever is larger, and similarly for the maximum.
}
\value{
  The robust range.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Werner Stahel, ETH Zurich}
\seealso{ \code{\link{plcoord}} }
\examples{
  x <- c(rnorm(20),rnorm(3,5,20))
  robrange(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}

