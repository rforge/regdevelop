\name{simresiduals}
\alias{simresiduals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate Residuals}
\description{
  Simulates residuals for a given normal regression model
}
\usage{
simresiduals(object, nrep, resgen=NULL, glm.restype="deviance")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{result of fitting a regression}
  \item{nrep}{number of replicates}
  \item{resgen}{if a function, it is used to generate random residuals,
    which will be multiplied by object\$sigma. If TRUE, \code{rnorm} will
    be used. If NULL (default) the standardized residuals of
    \code{object} will be randomly permuted}
  \item{glm.restype}{type of residuals to be generated (for glm)}
}
\details{
  The simulated residuals are obtained by replacing the response
  variable by permuted standardized residuals of the fitted regression,
  multiplied by the scale \code{object\$sigma}.
  If standarized residuals are not available, ordinary residuals are
  used. 
}
\value{
  A list consisting of
  \item{simres}{a matrix of which each column contains an set of
    simulated residuals}
  \item{simstres}{corresponding standardized residuals}
}
\author{Werner A. Stahel, ETH Zurich}
\examples{
"will follow"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
