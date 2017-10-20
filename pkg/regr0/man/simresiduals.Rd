\name{simresiduals}
\alias{simresiduals}
\alias{simresiduals.default}
\alias{simresiduals.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate Residuals}
\description{
  Simulates residuals for a given normal regression model
}
\usage{
simresiduals(object, ...)
\method{simresiduals}{default}(object, nrep, resgen=NULL)
\method{simresiduals}{glm}(object, nrep, resgen=NULL, glm.restype="working")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{result of fitting a regression}
  \item{nrep}{number of replicates}
  \item{resgen}{if a function, it is used to generate random residuals,
    which will be multiplied by object\$sigma. If TRUE, \code{rnorm} will
    be used. If NULL (default) the standardized residuals of
    \code{object} will be randomly permuted}
  \item{glm.restype}{type of residuals to be generated (for glm)
    Warning: type "deviance" may produce NAs.}
  \item{...}{further arguments passed to forthcoming methods.}
}
\details{
  The simulated residuals are obtained for the default method
  by replacing the response variable by permuted standardized residuals
  of the fitted regression, multiplied by the scale
  \code{object\$sigma}, then fitting the model to these residuals and
  getting the reseiduals  from this new fit.
  This is repeated \code{nrep} times.
  If standarized residuals are not available, ordinary residuals are
  used.

  For the \code{glm} method, the values of the response variable are
  obtained from simulating according to the model (binomial or Poisson),
  and the model is re-fitted to these generated values.
}
\value{
  A matrix of which each column contains an set of simulated residuals.
  If standardized residuals are available,
  attribute code{"stres"} is the matrix containing corresponding
  standardized residuals.
}
\author{Werner A. Stahel, ETH Zurich}
\examples{
"will follow"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
