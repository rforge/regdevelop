\name{polrregr}
\alias{polrregr}
\title{Fit a polr Regression and Return a Corresponding regr object
}
\description{
  Fit a regression for an ordered target variable by \code{polr}
  in a way adapted to the \code{regr} function.
}
\usage{
polrregr(formula, data, weights, start, ...,
         subset, na.action, contrasts = NULL, Hess = FALSE, model = TRUE,
         method = c("logistic", "probit", "cloglog", "cauchit"))
}
\arguments{
  \item{formula, data, weights, start, subset, na.action, contrasts,
    Hess, model, method, ...}{see \code{?\link[MASS]{polr}} of package \pkg{MASS}.
  }
}
\value{
  see \code{?\link[MASS]{polr}} of package \pkg{MASS} and \code{?\link{regr}}.
}
\author{Werner A. Stahel, ETH Zurich}

\seealso{
  \code{\link{polr}}, \code{\link{regr}}
}
\examples{
data(housing, package="MASS")
(t.r <- polrregr(Sat ~ Infl + Type + Cont, weights = housing$Freq,
                data = housing))
}
\keyword{regression}
