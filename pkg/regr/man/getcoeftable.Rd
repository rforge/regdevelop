\name{getcoeftable}
\alias{getcoeftable}
\alias{getscalepar}
\alias{getcoeffactor}

\title{Extract Componente of a Fit
}
\description{Retrieve the table of coefficients and standard errors,
  or the scale parameter, or the factors needed for standardizing
  coefficients from diverse model fitting results
}

\usage{
getcoeftable(object)
getscalepar(object)
getcoeffactor(object)
}

\arguments{
  \item{object}{results of a model fitting function
  }
}

\details{Object \code{regrModelClasses} contains the names of the
  classes for which the result should work.
  For other model classes, the function is not tested and may fail.
}

\value{For \code{getcoeftable}:
  Matrix containing at least the two columns containing the estimated
  coefficients (first column) and the standard errors (second column).\cr 
  For \code{getscalepar}: scale parameter. \cr
  For \code{getcoeffactor}: vector of factors, with attributes
  \code{scale}, \code{fitclass} and \code{family} or \code{dist}
  according to \code{object}.
}
\author{Werner A. Stahel}

\examples{
  rr <- lm(Fertility ~ . , data = swiss)
  getcoeftable(rr)
  getscalepar(rr)

  data(ovarian, package="survival")
  rs <- survival::survreg(survival::Surv(futime, fustat) ~ ecog.ps + rx,
                data = ovarian, dist = "weibull")
  getcoeftable(rs)
  getcoeffactor(rs)
}
\keyword{utilities}% use one of  RShowDoc("KEYWORDS")

