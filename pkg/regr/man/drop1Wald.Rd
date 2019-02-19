\name{drop1Wald}
\alias{drop1Wald}

\title{
Drop Single Terms of a Model and Calculate Respective Wald Tests
}
\description{
  'drop1Wald' calculates tests for single term deletions based on the
  covariance matrix of estimated coefficients instead of re-fitting a
  reduced model. This helps in cases where re-fitting is not feasible,
  inappropriate or costly.
}
\usage{
drop1Wald(object, scope=NULL, scale = 0, test = c("none", "Chisq", "F"),
  k = 2, ...)
}

\arguments{
  \item{object}{a fitted model of class 'regr'.
}
\item{scope}{a formula giving the terms to be considered for dropping.
  If 'NULL', 'drop.scope(object)' is obtained
}
  \item{scale}{an estimate of the residual mean square to be used in
          computing Cp. Ignored if '0' or 'NULL'.
}
  \item{test}{see \code{\link{drop1}}
}
  \item{k}{the penalty constant in AIC / Cp.
}
  \item{...}{further arguments, ignored}
}
\details{
  The test statistics and Cp and AIC values are calculated on the basis
  of the estimated coefficients and their (unscaled) covariance matrix
  as provided by the fitting object.
  The function may be used for all model fitting objects that contain
  these two components as \code{$coefficients} and \code{$cov.unscaled}.
}
\value{
   An object of class 'anova' summarizing the differences in fit
   between the models.
}
\author{
   Werner A. Stahel, Seminar for Statistics, ETH Zurich
}
\note{
  drop1Wald is used in \code{\link{regr}} for models of class 'lm' or
  'lmrob' for preparing the 'testcoef' table.
}

\seealso{
  \code{\link{drop1}}
}
\examples{
  data(d.blast)
  r.blast <-
       regr(log10(tremor)~location+log10(distance)+log10(charge),
            data=d.blast)
  drop1(r.blast)
  rr.blast <-
       regr(log10(tremor)~location+log10(distance)+log10(charge),
            data=d.blast, robust=TRUE)
  drop1(rr.blast)
}
\keyword{regression}

