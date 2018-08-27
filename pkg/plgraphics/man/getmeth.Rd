\name{getmeth}
\Rdversion{1.1}
\alias{getmeth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get S3 method of a generic function
}
\description{identical to \code{getS3method}}
\usage{
getmeth(fn, mt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{name of generic function, quoted or unquoted}
  \item{mt}{name of method, quoted or unquoted}
}
\value{Source code of the method
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{getS3method}}}
\examples{getmeth(simresiduals, glm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{debugging}
