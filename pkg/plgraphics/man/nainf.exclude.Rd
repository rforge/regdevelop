\name{nainf.exclude}
\alias{nainf.exclude}
\alias{nainf.omit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Drop Rows Containing NA or Inf}
\description{
  Drops the rows of a data frame that contain an NA, an NaN, or an Inf value
}
\usage{
nainf.exclude(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an R object, typically a data frame}
  \item{\dots}{further arguments special methods could require.}
}
\details{
  This is a simple modification of \code{\link{na.omit}} and
    \code{\link{na.exclude}}
}
\value{
  The value is of the same type as the argument \code{object},
  with possibly less elements.
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{na.omit}}}
\examples{
t.d <- data.frame(V1=c(1,2,NA,4), V2=c(11,12,13,Inf))
nainf.exclude(t.d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{NA}
\keyword{manip}% __ONLY ONE__ keyword per line
