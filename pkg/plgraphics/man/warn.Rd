\name{warn}
\alias{warn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{List Warnings}
\description{
  Gives a List of Warnings
}
\usage{
warn()
}
%- maybe also 'usage' for other objects documented here.
\details{
  This function simplyfies the output of \code{\link{warnings}} if there
  are several identical warnings, by counting their occurence
}
\value{
  the table of warnings
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{warnings}}}
\examples{
##  not to be run as an example()
##  for (i in 3:6)
##    m <- matrix(1:7, 3,i)
##  warn()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
