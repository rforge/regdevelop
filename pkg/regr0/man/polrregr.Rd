\name{polrregr}
\Rdversion{1.1}
\alias{polrregr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a polr Regression and Return a Corresponding regr object
}
\description{
  Fit a regression for an ordered target variable by \code{polr}
  and return a corresponding regr object.
}
\usage{
polrregr(formula, data, weights, start, ..., 
subset, na.action, contrasts = NULL, Hess = FALSE, model = TRUE,
method = c("logistic", "probit", "cloglog", "cauchit"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula, data, weights, start, subset, na.action, contrasts,
    Hess, model, method, ...}{see \code{?polr} of package MASS
  }
}
\value{
  see \code{?polr} of package MASS and \code{?regr}
}
\author{Werner A. Stahel, ETH Zurich}

\seealso{
\code{\link{polr}}, \code{\link{regr}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
